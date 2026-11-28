"""Generate the synthetic long-range test sets (deterministic).

  silesia_snap.bin : 4 "backup snapshots" of the Silesia tar; each snapshot
                     is the previous one with ~3000 random local edits
                     (overwrite / insert / delete, 1-256 bytes).
  db_snap.sql      : 3 daily dumps of a synthetic shop database in mysqldump
                     style; between dumps ~3% of rows change and ~4% are added.

usage: python gen_sets.py SILESIA_DIR OUTDIR
"""
import os
import random
import sys


def silesia_snapshots(src_dir, out, nsnap=4, edits=3000, seed=1):
    names = sorted(os.listdir(src_dir))
    base = bytearray()
    for n in names:
        base += open(os.path.join(src_dir, n), "rb").read()
    rnd = random.Random(seed)
    snap = bytes(base)
    with open(out, "wb") as f:
        for k in range(nsnap):
            f.write(snap)
            # apply edits back to front so positions stay valid
            pos = sorted(rnd.randrange(len(snap)) for _ in range(edits))
            parts, last = [], 0
            for p in pos:
                if p < last:
                    continue
                parts.append(snap[last:p])
                op = rnd.random()
                ln = rnd.randint(1, 256)
                if op < 0.4:      # overwrite
                    parts.append(rnd.randbytes(ln)); last = p + ln
                elif op < 0.7:    # insert (copy of other content)
                    q = rnd.randrange(len(snap) - ln)
                    parts.append(snap[q:q + ln]); last = p
                else:             # delete
                    last = p + ln
            parts.append(snap[last:])
            snap = b"".join(parts)


def db_snapshots(vocab_file, out, nsnap=3, seed=2):
    rnd = random.Random(seed)
    words = sorted({w.strip(".,;:!?\"'()[]").lower() for w in open(vocab_file, "rb").read().decode("latin-1").split()
                    if w.isalpha() and 3 <= len(w) <= 12})
    first = [w.capitalize() for w in words[:4000]]
    last = [w.capitalize() for w in words[4000:12000]]
    cities = [w.capitalize() for w in words[12000:12600]]
    status = ["active", "inactive", "pending", "banned"]
    ostat = ["new", "paid", "shipped", "delivered", "returned", "cancelled"]

    def user(i):
        fn, ln = rnd.choice(first), rnd.choice(last)
        return [i, fn, ln, f"{fn.lower()}.{ln.lower()}{rnd.randint(1, 999)}@example.{rnd.choice(['com', 'net', 'org'])}",
                f"2{rnd.randint(18, 25):03d}-{rnd.randint(1, 12):02d}-{rnd.randint(1, 28):02d} {rnd.randint(0, 23):02d}:{rnd.randint(0, 59):02d}:{rnd.randint(0, 59):02d}",
                f"{rnd.randint(0, 99999)}.{rnd.randint(0, 99):02d}", rnd.choice(cities), rnd.choice(status)]

    def order(i, nusers):
        return [i, rnd.randrange(1, nusers), rnd.randint(1, 5000), rnd.randint(1, 9),
                f"{rnd.randint(1, 2000)}.{rnd.randint(0, 99):02d}", rnd.choice(ostat),
                f"2{rnd.randint(18, 25):03d}-{rnd.randint(1, 12):02d}-{rnd.randint(1, 28):02d}",
                " ".join(rnd.choice(words) for _ in range(rnd.randint(0, 6)))]

    users = {i: user(i) for i in range(1, 600001)}
    orders = {i: order(i, 600000) for i in range(1, 1500001)}

    def fmt(v):
        return str(v) if isinstance(v, int) else "'" + v.replace("'", "''") + "'"

    def dump(f, day):
        f.write(f"-- MySQL dump, snapshot day {day}\n--\n".encode())
        for name, tab, cols in (("users", users, "id,first,last,email,created,balance,city,status"),
                                ("orders", orders, "id,user_id,product_id,qty,amount,status,date,note")):
            f.write(f"DROP TABLE IF EXISTS `{name}`;\nCREATE TABLE `{name}` ({cols});\n".encode())
            keys = sorted(tab)
            for s in range(0, len(keys), 500):
                rows = ",".join("(" + ",".join(fmt(v) for v in tab[k]) + ")" for k in keys[s:s + 500])
                f.write(f"INSERT INTO `{name}` VALUES {rows};\n".encode())

    with open(out, "wb") as f:
        for day in range(nsnap):
            dump(f, day)
            nu, no = max(users), max(orders)
            for _ in range(len(users) * 3 // 100):
                k = rnd.randrange(1, nu)
                if k in users:
                    users[k][5] = f"{rnd.randint(0, 99999)}.{rnd.randint(0, 99):02d}"; users[k][7] = rnd.choice(status)
            for _ in range(len(orders) * 3 // 100):
                k = rnd.randrange(1, no)
                if k in orders:
                    orders[k][5] = rnd.choice(ostat)
            for _ in range(len(orders) // 200):
                orders.pop(rnd.randrange(1, no), None)
            for i in range(nu + 1, nu + 1 + len(users) * 4 // 100):
                users[i] = user(i)
            for i in range(no + 1, no + 1 + len(orders) * 4 // 100):
                orders[i] = order(i, nu)


if __name__ == "__main__":
    sdir, odir = sys.argv[1], sys.argv[2]
    silesia_snapshots(sdir, os.path.join(odir, "silesia_snap.bin"))
    db_snapshots(os.path.join(sdir, "dickens"), os.path.join(odir, "db_snap.sql"))
