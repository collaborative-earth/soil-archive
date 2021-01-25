import pandas as pd
from typing import List


def combine_WISE_site_data(
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


def compute_WISE_carbon_content(df: pd.DataFrame, depth: float = 100.0):
    """Given a DataFrame, which is assumed to be the output of combine_site_data,
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


def _filter_ISCN(df: pd.DataFrame) -> pd.DataFrame:
    """Given a raw ISCN layer-level DataFrame, clean variable names and filter
    to only layers from valid profiles."""
    filtered_df = df[
        df["soc (g cm-2)"].notnull()
        & df["layer_top (cm)"].notnull()
        & df["layer_bot (cm)"].notnull()
    ][
        [
            "dataset_name_sub",
            "dataset_name_soc",
            "lat (dec. deg)",
            "long (dec. deg)",
            "datum (datum)",
            "state (state_province)",
            "country (country)",
            "site_name",
            "observation_date (YYYY-MM-DD)",
            "profile_name",
            "layer_name",
            "layer_top (cm)",
            "layer_bot (cm)",
            "soc (g cm-2)",
            "soc_carbon_flag",
        ]
    ].reset_index(
        drop=True
    )
    filtered_df.rename(
        {
            "lat (dec. deg)": "latitude",
            "long (dec. deg)": "longitude",
            "datum (datum)": "datum",
            "state (state_province)": "state",
            "country (country)": "country",
            "observation_date (YYYY-MM-DD)": "observation_date",
            "layer_top (cm)": "layer_top",
            "layer_bot (cm)": "layer_bot",
            "soc (g cm-2)": "soc",
        },
        axis=1,
        inplace=True,
    )

    return filtered_df


def _verify_ISCN_quality(df) -> [pd.DataFrame, pd.DataFrame]:
    grouped_df = (
        df.groupby(
            by=[
                "dataset_name_sub",
                "dataset_name_soc",
                "latitude",
                "longitude",
                "site_name",
                "observation_date",
                "profile_name",
            ]
        )
        .agg(
            {
                "scaled_soc": "sum",
                "full_length": "sum",
                "truncated_length": "sum",
                "layer_top": "min",
                "layer_bot": "max",
                "layer_bot_trunc": "max",
            }
        )
        .reset_index()
    )

    msk = (
        grouped_df.truncated_length <= grouped_df.layer_bot_trunc - grouped_df.layer_top
    )
    grouped_df = grouped_df[
        [
            "dataset_name_sub",
            "dataset_name_soc",
            "latitude",
            "longitude",
            "site_name",
            "observation_date",
            "profile_name",
            "scaled_soc",
            "truncated_length",
            "layer_top",
            "layer_bot",
            "layer_bot_trunc",
        ]
    ].rename(
        {
            "layer_top": "min_depth",
            "layer_bot": "max_depth",
            "layer_bot_trunc": "max_depth_trunc",
        },
        axis="columns",
    )
    return grouped_df[msk], grouped_df[~msk]


def compute_ISCN_carbon_content(df: pd.DataFrame, depth: float = 100.0):
    df = _filter_ISCN(df)

    # Truncate any segment which extends below user-specified depth
    df["layer_bot_trunc"] = df.layer_bot.clip(upper=depth)
    df["truncated_length"] = df.layer_bot_trunc - df.layer_top
    df["full_length"] = df.layer_bot - df.layer_top

    # Filter to segments which have non-zero truncated length and scale computed
    # carbon content of layer by the ratio of the truncated length to full length
    # This step assumes distribution of carbon over layer is ~constant.
    df = df[df.truncated_length.gt(0)]
    df["scaled_soc"] = (df.truncated_length / df.full_length) * df.soc * 10 ** 4

    valid_df, invalid_df = _verify_dataset_quality(df)

    return (
        valid_df[
            [
                "dataset_name_sub",
                "dataset_name_soc",
                "latitude",
                "longitude",
                "site_name",
                "observation_date",
                "profile_name",
                "scaled_soc",
                "truncated_length",
                "min_depth",
            ]
        ]
        .rename(
            {"scaled_soc": "carbon_stock", "truncated_length": "depth"}, axis="columns"
        )
        .reset_index(drop=True),
        invalid_df.reset_index(drop=True),
    )
