import pandas as pd

from src.splits.leakage_safe_split import build_strict_forward_site_grouped_splits


def test_build_strict_forward_site_grouped_splits_has_no_subject_overlap():
    df = pd.DataFrame(
        {
            "subject_id": ["s1", "s2", "s3", "s4", "s5", "s6"],
            "site_id": ["A", "A", "B", "B", "C", "C"],
            "collection_time": pd.to_datetime(
                [
                    "2020-01-01",
                    "2020-02-01",
                    "2020-03-01",
                    "2020-04-01",
                    "2020-05-01",
                    "2020-06-01",
                ]
            ),
        }
    )

    splits = build_strict_forward_site_grouped_splits(
        df,
        time_col="collection_time",
        subject_col="subject_id",
        train_frac=0.5,
        valid_frac=0.25,
    )

    train_subjects = set(df.loc[splits["train_idx"], "subject_id"])
    valid_subjects = set(df.loc[splits["valid_idx"], "subject_id"])
    test_subjects = set(df.loc[splits["test_idx"], "subject_id"])

    assert train_subjects.isdisjoint(valid_subjects)
    assert train_subjects.isdisjoint(test_subjects)
    assert valid_subjects.isdisjoint(test_subjects)
    assert len(splits["train_idx"]) > 0
    assert len(splits["valid_idx"]) > 0
    assert len(splits["test_idx"]) > 0
