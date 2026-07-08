"""
Download Level 2 Sentinel-3 products from the Copernicus Data Space Ecosystem
(CDSE) for a given lon/lat point and date range.

Meant to be imported (not run standalone):

    from download_sentinel3 import download_sentinel3_l2

    files = download_sentinel3_l2(
        lon=-4.5, lat=48.4,
        start_date="2024-06-01", end_date="2024-06-30",
        out_dir="./data/sentinel3",
        product_type="OL_2_WFR___",   # OLCI full-resolution water product
    )

Credentials: create a free account at https://dataspace.copernicus.eu and
either export COPERNICUS_USERNAME / COPERNICUS_PASSWORD as environment
variables, or pass username=/password= explicitly. Never hardcode
credentials in a script that gets committed to git.

Common Sentinel-3 Level 2 product types (used for the `product_type` /
`contains(Name, ...)` filter):
    OL_2_WFR___   OLCI full resolution, water (ocean colour)      <- default
    OL_2_WRR___   OLCI reduced resolution, water
    SL_2_WST___   SLSTR sea surface temperature
    SL_2_LST___   SLSTR land surface temperature
    SR_2_WAT___   SRAL altimetry, water
"""

from __future__ import annotations

import os
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

import requests

IDENTITY_URL = (
    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/"
    "protocol/openid-connect/token"
)
CATALOGUE_URL = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
DOWNLOAD_URL = "https://zipper.dataspace.copernicus.eu/odata/v1/Products"


def get_access_token(username: str, password: str) -> dict:
    """Exchange CDSE credentials for an access/refresh token pair."""
    response = requests.post(
        IDENTITY_URL,
        data={
            "client_id": "cdse-public",
            "grant_type": "password",
            "username": username,
            "password": password,
        },
        timeout=30,
    )
    if not response.ok:
        raise RuntimeError(
            f"CDSE authentication failed ({response.status_code}): {response.text}"
        )
    return response.json()


def refresh_access_token(refresh_token: str) -> dict:
    """Use a refresh token to obtain a new access token without re-sending credentials."""
    response = requests.post(
        IDENTITY_URL,
        data={
            "client_id": "cdse-public",
            "grant_type": "refresh_token",
            "refresh_token": refresh_token,
        },
        timeout=30,
    )
    if not response.ok:
        raise RuntimeError(
            f"CDSE token refresh failed ({response.status_code}): {response.text}"
        )
    return response.json()


def _to_iso(date_str: str, end_of_day: bool = False) -> str:
    """Convert a 'YYYY-MM-DD' (or full ISO) string to an OData-friendly UTC timestamp."""
    if "T" in date_str:
        return date_str if date_str.endswith("Z") else date_str + "Z"
    suffix = "T23:59:59.999Z" if end_of_day else "T00:00:00.000Z"
    # validates the date is well-formed before we build the filter string
    datetime.strptime(date_str, "%Y-%m-%d")
    return f"{date_str}{suffix}"


def search_sentinel3_l2(
    lon: float,
    lat: float,
    start_date: str,
    end_date: str,
    product_type: str = "OL_2_WFR___",
    max_records: int = 100,
) -> list[dict]:
    """
    Query the CDSE OData catalogue for Sentinel-3 Level 2 products whose
    footprint contains (lon, lat) and whose sensing time falls in
    [start_date, end_date].

    Returns a list of dicts with keys: id, name, size, sensing_start.
    """
    start_iso = _to_iso(start_date)
    end_iso = _to_iso(end_date, end_of_day=True)

    filter_str = (
        f"Collection/Name eq 'SENTINEL-3' "
        f"and contains(Name,'{product_type}') "
        f"and OData.CSC.Intersects(area=geography'SRID=4326;POINT({lon} {lat})') "
        f"and ContentDate/Start gt {start_iso} "
        f"and ContentDate/Start lt {end_iso}"
    )

    products: list[dict] = []
    skip = 0
    page_size = min(max_records, 1000)
    while True:
        response = requests.get(
            CATALOGUE_URL,
            params={
                "$filter": filter_str,
                "$orderby": "ContentDate/Start asc",
                "$top": page_size,
                "$skip": skip,
            },
            timeout=60,
        )
        if not response.ok:
            raise RuntimeError(
                f"CDSE search failed ({response.status_code}): {response.text}"
            )
        page = response.json().get("value", [])
        for item in page:
            products.append(
                {
                    "id": item["Id"],
                    "name": item["Name"],
                    "size": item.get("ContentLength"),
                    "sensing_start": item.get("ContentDate", {}).get("Start"),
                }
            )
        if len(page) < page_size or len(products) >= max_records:
            break
        skip += page_size

    return products[:max_records]


def download_product(
    product_id: str,
    product_name: str,
    out_dir: str | Path,
    access_token: str,
    refresh_token: Optional[str] = None,
    chunk_size: int = 1024 * 1024,
) -> Path:
    """
    Download a single product by its CDSE Id to `out_dir/<product_name>.zip`,
    refreshing the access token if it expires mid-download.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{product_name}.zip"

    if out_path.exists():
        return out_path

    url = f"{DOWNLOAD_URL}({product_id})/$value"
    headers = {"Authorization": f"Bearer {access_token}"}

    with requests.Session() as session:
        response = session.get(url, headers=headers, stream=True, timeout=60)

        if response.status_code == 401 and refresh_token:
            token_data = refresh_access_token(refresh_token)
            headers["Authorization"] = f"Bearer {token_data['access_token']}"
            response = session.get(url, headers=headers, stream=True, timeout=60)

        if not response.ok:
            raise RuntimeError(
                f"Download failed for {product_name} "
                f"({response.status_code}): {response.text[:500]}"
            )

        tmp_path = out_path.with_suffix(".part")
        with open(tmp_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
        tmp_path.rename(out_path)

    return out_path


def download_sentinel3_l2(
    lon: float,
    lat: float,
    start_date: str,
    end_date: str,
    out_dir: str | Path,
    product_type: str = "OL_2_WFR___",
    username: Optional[str] = None,
    password: Optional[str] = None,
    max_records: int = 100,
) -> list[Path]:
    """
    Search and download all Sentinel-3 Level 2 products of `product_type`
    covering (lon, lat) within [start_date, end_date].

    Credentials fall back to the COPERNICUS_USERNAME / COPERNICUS_PASSWORD
    environment variables if `username`/`password` are not given.

    Dates may be 'YYYY-MM-DD' or full ISO timestamps.

    Returns the list of local file paths downloaded (or already present).
    """
    username = username or os.environ.get("COPERNICUS_USERNAME")
    password = password or os.environ.get("COPERNICUS_PASSWORD")
    if not username or not password:
        raise ValueError(
            "Copernicus credentials not found. Pass username=/password= or "
            "set the COPERNICUS_USERNAME / COPERNICUS_PASSWORD environment "
            "variables (sign up at https://dataspace.copernicus.eu)."
        )

    products = search_sentinel3_l2(
        lon, lat, start_date, end_date,
        product_type=product_type, max_records=max_records,
    )
    if not products:
        print(f"No {product_type} products found for ({lon}, {lat}) "
              f"between {start_date} and {end_date}.")
        return []

    token_data = get_access_token(username, password)
    access_token = token_data["access_token"]
    refresh_token = token_data.get("refresh_token")
    token_fetched_at = time.time()

    downloaded = []
    for product in products:
        # CDSE access tokens expire after ~10 minutes; refresh proactively
        # if a long search/download loop has been running for a while.
        if time.time() - token_fetched_at > 540 and refresh_token:
            token_data = refresh_access_token(refresh_token)
            access_token = token_data["access_token"]
            refresh_token = token_data.get("refresh_token", refresh_token)
            token_fetched_at = time.time()

        print(f"Downloading {product['name']} ...")
        path = download_product(
            product["id"], product["name"], out_dir,
            access_token, refresh_token=refresh_token,
        )
        downloaded.append(path)

    return downloaded
