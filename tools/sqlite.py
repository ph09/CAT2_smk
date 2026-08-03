"""
Tools to wrap around the native sqlite package. Necessary due to bugs in how pandas interacts with sqlalchemy.
"""
import fcntl
import sqlite3 as sql
import time
import warnings

__author__ = "Ian Fiddes"


class ExclusiveSqlConnection(object):
    """Context manager for an exclusive SQL connection"""
    def __init__(self, path, timeout=6000):
        self.path = path
        self.timeout = timeout

    def __enter__(self):
        self.con = sql.connect(self.path, timeout=self.timeout, isolation_level="EXCLUSIVE")
        try:
            self.con.execute("BEGIN EXCLUSIVE")
        except sql.OperationalError:
            raise RuntimeError("Database still locked after {} seconds.".format(self.timeout))
        return self.con

    def __exit__(self, exception_type, exception_val, trace):
        try:
            if exception_type is None:
                self.con.commit()
            else:
                self.con.rollback()
        finally:
            self.con.close()


class _DeferredCommitConnection(object):
    """Proxy around sqlite3.Connection that ignores commit()/close() from pandas.

    pandas.io.sql.SQLiteDatabase.run_transaction() always calls con.commit()
    after writing one table. That releases SQLite's exclusive lock and lets a
    concurrent evaluate job interleave mid-replace (DROP/CREATE without INSERT),
    which leaves empty or missing metrics/evaluation tables.

    On Python 3.12+, Connection.commit is read-only, so we cannot monkeypatch it
    on the real connection; this proxy is required.
    """

    def __init__(self, con):
        object.__setattr__(self, "_con", con)

    def commit(self):
        return None

    def close(self):
        # Owner closes the real connection.
        return None

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, "_con"), name)


def write_dataframes(path, table_frames, timeout=6000, retries=12, retry_sleep=2.0):
    """Write many DataFrames into one SQLite DB under a single exclusive transaction.

    Holds a sidecar flock for the whole batch, begins BEGIN EXCLUSIVE once, and
    defers pandas' per-table commits so concurrent writers cannot observe a
    half-replaced table.
    """
    path = str(path)
    lock_path = path + ".write.lock"
    last_err = None

    for attempt in range(retries):
        try:
            with open(lock_path, "w") as lockf:
                fcntl.flock(lockf, fcntl.LOCK_EX)
                # isolation_level=None: we manage the transaction explicitly.
                con = sql.connect(path, timeout=timeout, isolation_level=None)
                try:
                    con.execute("BEGIN EXCLUSIVE")
                    wrapped = _DeferredCommitConnection(con)
                    with warnings.catch_warnings():
                        # Proxy is not a raw sqlite3.Connection; pandas warns but
                        # still uses the DBAPI2 path correctly via __getattr__.
                        warnings.simplefilter("ignore", UserWarning)
                        for table_name, df in table_frames:
                            df.to_sql(table_name, wrapped, if_exists="replace", index=True)
                    con.commit()
                except Exception:
                    try:
                        con.rollback()
                    except sql.Error:
                        pass
                    raise
                finally:
                    con.close()
            return
        except sql.OperationalError as e:
            last_err = e
            if "locked" not in str(e).lower() or attempt + 1 >= retries:
                raise
            time.sleep(retry_sleep * (attempt + 1))

    raise last_err


def attach_database(con, path, name):
    """
    Attaches another database found at path to the name given in the given connection.
    """
    con.execute("ATTACH DATABASE '{}' AS {}".format(path, name))


def open_database(path, timeout=6000):
    """opens a database, returning the connection and cursor objects."""
    con = sql.connect(path, timeout=timeout)
    cur = con.cursor()
    return con, cur
