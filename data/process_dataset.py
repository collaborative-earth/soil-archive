import pandas as pd
from typing import List


def _combine_site_data(
    horizon_path: str,
    site_path: str,
    cols: List[str] = [
        "WISE3_ID",  # Unique ID for soil in question (ISO code + number)
        "HONU",  #  sequential horizon number
        "DESIG",  # horizon designation code
        "TOPDEP",  # upper depth of horizon (cm)
        "BOTDEP",  # bottom depth of horizon (cm)
        "ORGC",  # organic carbon content (g kg^-1)
        "BULKDENS",  # bulk density (g cm^-3)
        "LONDD",  # longitude
        "LATDD",  # latitude
        "LONLAT_ACC",  # accuracy of location
    ],
) -> pd.DataFrame:
    """Given paths to horizon and site data, join the two and return a DataFrame
    containing the specified columns."""
    horizon_df = pd.read_csv(horizon_path)
    site_df = pd.read_csv(site_path)
    combined_df = horizon_df.merge(
        site_df, how="inner", left_on="WISE3_ID", right_on="WISE3_id"
    )

    return combined_df[cols]


def _compute_carbon_content(df: pd.DataFrame, depth: float = 100.0):
    """Given a DataFrame, which is assumed to be the output of _combine_site_data,
    compute the carbon stock up to 1m depth. Places where data is not available
    up to 1m will have a depth value reported which is less than 100 cm."""
    df = df.copy()
    df["BOTDEP_trunc"] = df["BOTDEP"].clip(upper=depth)
    df["segment_length"] = df["BOTDEP_trunc"] - df["TOPDEP"]

    df = df[df.BULKDENS.notnull() & df.segment_length.gt(0)]
    df["carbon_stock"] = 10 * df["ORGC"] * df["BULKDENS"] * df["segment_length"]

    return (
        df.groupby(by=["WISE3_ID", "LONDD", "LATDD", "LONLAT_ACC"])
        .agg({"carbon_stock": "sum", "segment_length": "sum"})
        .rename({"segment_length": "depth"}, axis="columns")
        .reset_index()
    )
