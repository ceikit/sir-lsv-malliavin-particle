"""
Pull one day of CME EUR/USD options-on-futures settlements from Databento and
write them to a single CSV (the input for the calibration test).

What it requests (dataset GLBX.MDP3, parent symbology):
  * definitions  for EUU.OPT (monthly/quarterly European EUR/USD options) and 6E.FUT
  * statistics   for the same symbols: settlement price, open interest, cleared volume

Before downloading it prints Databento's own cost estimate and asks for confirmation.

Setup (once):   pip install databento
                export DATABENTO_API_KEY=db-...        (do not paste the key into a chat)
Run:            python data/pull_eurusd_options_databento.py 2026-11-12
Output:         eurusd_cme_options_2026-11-12.csv   (licensed data: do not commit it)

Not yet run against the live API: the request layout follows Databento's documented
examples and the installed client's signatures, but the symbol root EUU and the
date conventions are unverified.

Known issue: build_table keeps the last statistics message per instrument and does not
look at update_action, so a statistic that was published and then deleted is still used.
"""
import sys
import pandas as pd

DATASET = "GLBX.MDP3"
SYMBOLS = ["EUU.OPT", "6E.FUT"]


def build_table(defs: pd.DataFrame, stats: pd.DataFrame, trade_date, stat_types) -> pd.DataFrame:
    """Join the last definition of each instrument with its statistics for trade_date."""
    defs = defs.sort_index().drop_duplicates("instrument_id", keep="last")
    defs = defs[defs["instrument_class"].astype(str).isin(["C", "P", "F"])]       # drop spreads
    keep = [c for c in ["instrument_id", "raw_symbol", "instrument_class", "strike_price",
                        "expiration", "underlying", "asset", "unit_of_measure_qty"] if c in defs.columns]
    out = defs[keep].copy()

    stats = stats.sort_index()
    stats = stats[pd.to_datetime(stats["ts_ref"]).dt.date == trade_date]
    for name, (stype, col) in stat_types.items():
        s = stats[stats["stat_type"] == stype].drop_duplicates("instrument_id", keep="last")   # last message wins
        cols = {col: name}
        if name == "settlement":
            cols["stat_flags"] = "settlement_flags"
        out = out.merge(s[["instrument_id"] + list(cols)].rename(columns=cols), on="instrument_id", how="left")

    out.insert(0, "trade_date", str(trade_date))
    out = out[out["settlement"].notna()]
    return out.sort_values(["instrument_class", "expiration", "strike_price"]).reset_index(drop=True)


def main():
    import databento as db

    date = pd.Timestamp(sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else "2026-11-12")
    start, end_def, end_stat = date, date + pd.Timedelta(days=1), date + pd.Timedelta(days=2)
    client = db.Historical()                          # reads DATABENTO_API_KEY

    req = dict(dataset=DATASET, symbols=SYMBOLS, stype_in="parent")
    cost = (client.metadata.get_cost(schema="definition", start=start, end=end_def, **req)
            + client.metadata.get_cost(schema="statistics", start=start, end=end_stat, **req))
    print(f"Databento cost estimate for this pull: ${cost:.4f}")
    if "--yes" not in sys.argv and input("Download? [y/N] ").strip().lower() != "y":
        sys.exit("cancelled")

    defs = client.timeseries.get_range(schema="definition", start=start, end=end_def, **req).to_df()
    stats = client.timeseries.get_range(schema="statistics", start=start, end=end_stat, **req).to_df()
    print(f"definitions: {len(defs)} rows, statistics: {len(stats)} rows")

    table = build_table(defs, stats, date.date(), {
        "settlement": (db.StatType.SETTLEMENT_PRICE, "price"),
        "open_interest": (db.StatType.OPEN_INTEREST, "quantity"),
        "cleared_volume": (db.StatType.CLEARED_VOLUME, "quantity"),
    })
    path = f"eurusd_cme_options_{date.date()}.csv"
    table.to_csv(path, index=False)
    n_opt = int(table["instrument_class"].astype(str).isin(["C", "P"]).sum())
    print(f"wrote {path}: {n_opt} options, {len(table) - n_opt} futures")
    if n_opt:
        opt = table[table["instrument_class"].astype(str).isin(["C", "P"])]
        print("option expiries:", ", ".join(sorted({str(pd.Timestamp(e).date()) for e in opt["expiration"]})))


if __name__ == "__main__":
    main()
