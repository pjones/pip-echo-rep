#!/usr/bin/env python
import json
import os.path
import subprocess

def download_urls(json_file: str, directory: str):
    with open(json_file) as file:
      urls = json.load(file)

    for link in urls:
        path = os.path.join(directory, link["name"])
        if not os.path.exists(path):
            print(path)
            subprocess.run(["curl", "-o", path, link["href"]])

if __name__ == "__main__":
    download_urls("data/PXD057758/links.json", "data/PXD057758")
