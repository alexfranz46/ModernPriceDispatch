import requests
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.parse import quote

### NOTE: Taken from Copilot and not yet tested.


# =========================
# CONFIGURATION
# =========================

CONTAINER_URL = "https://emidatasets.blob.core.windows.net/publicdata"

PREFIX = (
    "Datasets/Wholesale/"
    "DispatchAndPricing/"
    "DispatchEnergyPrices/"
)

DOWNLOAD_FOLDER = Path("DispatchEnergyPrices")


# =========================
# GET ALL FILES
# =========================

def get_all_blob_names(prefix):
    """
    Retrieve all file paths from Azure Blob Storage.

    Azure returns a maximum of 5000 results per request.
    If there are more files, the XML response contains a
    NextMarker value which we use to fetch the next page.
    """

    all_files = []
    marker = None

    while True:

        url = (
            f"{CONTAINER_URL}"
            f"?restype=container"
            f"&comp=list"
            f"&prefix={quote(prefix)}"
        )

        if marker:
            url += f"&marker={quote(marker)}"

        print(f"Querying: {url}")

        response = requests.get(url)
        response.raise_for_status()

        root = ET.fromstring(response.text)

        # Extract file names
        for blob in root.findall(".//Blob"):
            name = blob.find("Name").text
            if name.endswith(".csv"):
                all_files.append(name)

        # Check if there is another page
        next_marker_element = root.find(".//NextMarker")

        if (
            next_marker_element is None
            or not next_marker_element.text
        ):
            break

        marker = next_marker_element.text

        print(
            f"Found {len(all_files)} files so far..."
        )

    return all_files


# =========================
# DOWNLOAD A SINGLE FILE
# =========================

def download_file(blob_name):
    """
    Download a single blob.
    """

    if not blob_name.endswith(".csv"):
        return

    local_path = DOWNLOAD_FOLDER / Path(blob_name).relative_to(PREFIX)

    local_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    if local_path.exists():
        print(f"Skipping existing file: {local_path.name}")
        return

    file_url = f"{CONTAINER_URL}/{blob_name}"

    print(f"Downloading: {local_path.name}")

    response = requests.get(
        file_url,
        stream=True
    )

    response.raise_for_status()

    with open(local_path, "wb") as f:
        for chunk in response.iter_content(
            chunk_size=1024 * 1024
        ):
            if chunk:
                f.write(chunk)


# =========================
# MAIN PROGRAM
# =========================

def main():

    DOWNLOAD_FOLDER.mkdir(
        exist_ok=True
    )

    print("Getting file list...")

    files = get_all_blob_names(PREFIX)

    print(
        f"\nFound {len(files)} files.\n"
    )

    # TEST RUN: Only keep first 3 files
    # files = files[:3]

    for i, file in enumerate(files, start=1):

        print(
            f"[{i}/{len(files)}]"
        )

        download_file(file)

    print("\nFinished.")


if __name__ == "__main__":
    main()