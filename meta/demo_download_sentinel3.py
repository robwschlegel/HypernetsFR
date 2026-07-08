"""
Minimal demo: source download_sentinel3.py and fetch Sentinel-3 Level 2
products for a single point.

Requires COPERNICUS_USERNAME / COPERNICUS_PASSWORD to be set in the
environment (see download_sentinel3.py for details).
"""

from download_sentinel3 import download_sentinel3_l2

files = download_sentinel3_l2(
    lon=-1.04,
    lat=45.55,
    start_date="2024-06-01",
    end_date="2024-06-30",
    out_dir="../data/S3",
    product_type="OL_2_WFR___",  # OLCI full-resolution water product
)

print(f"Downloaded {len(files)} file(s):")
for f in files:
    print(f" - {f}")
