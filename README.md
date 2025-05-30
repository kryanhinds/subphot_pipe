---
# subphot\_pipe
---
**subphot\_pipe** is an automated image subtraction and photometry pipeline based on PSF cross-convolution. It supports subtraction relative to reference images from **Pan-STARRS**, **SDSS**, or a custom input image. Designed for transient and variable object analysis, it provides streamlined image alignment, PSF matching, and photometric measurements with extensive configurability.


## 🔧 Installation

Clone the repository:

```bash
git clone https://github.com/kryanhinds/subphot_pipe.git
cd subphot_pipe
```

### Set Up Conda Environment

For **macOS**:

```bash
conda env create -f environment.yml
```

For **Linux**:

```bash
conda env create -f environment_linux.yml
```

Then activate the environment:

```bash
conda activate subphot_pipe
```

---

## ⚙️ Configuration

### 1. Credentials File

Rename and edit the credentials template:

```bash
cp subphot_credentials_blank.py subphot_credentials.py
```

Update:

* `data1_path`: Path to your data directory (input images and output files will be handled here).
* `token`: Your [Fritz](https://fritz.science/) SkyPortal API token for uploading/retrieving photometry.
* `smtp2go_api_key` and `outlook_account`: For email notifications upon completion.
* `swarp_path`, `sex_path`, `sexpath`, `panstamps_path`: Paths to the required software installations.


### 2. Config Files

Update all config files in `subphot_pipe/config/` to reflect your current working directory:

* `config_files/prepsfex.sex`
* `config_files/align_sex.config`
* `config_files/config.swarp`
* `config_files/config_comb.swarp`
* `config_files/config_resize.swarp`
* `config_files/psfex_conf.psfex`
* `config_files/sex.config`
* `config_files/swarp_sdss.conf`


---

## 🔄 Parameters of Interest

* `image_size`: Used during image alignment.
* `starscale`, `search_rad`: Threshold-based star detection.
* `store_lc_ims`: Enables organizing images by target for better tracking.
* `proposals_arc`: Dictionary for scraping Quicklook/RecentData images based on proposal IDs.

---

## 🚀 Usage

### Basic Commands

Run on a single image:

```bash
python3 subphot_subtract.py -i h_e_20220206_1_1_1.fits
```

Photometry, cutouts and any logs will be saved to data1_path/photometry, data1_path/photometry/cut_outs and data1_path/photometry/<name>.log, respecitvely.

Run on all images in a folder:

```bash
python3 subphot_subtract.py -f SNDataFolder
```

Run on specific images:

```bash
python3 subphot_subtract.py -i image1.fits image2.fits image3.fits
```

Run on a wildcard and stack:

```bash
python3 subphot_subtract.py -i h_e_20220206_*.fits -stk
```

Advanced example with full options:

```bash
python3 subphot_subtract.py -f SNDataFolder -sn SNXXXX -b r i -cut -o SNDataResults -fp RAX DECX -log SNXXXX
```
This will run the code on all images in the folder SNDataFolder in the data1_path folder, for SN SNXXXX, in filters r and i, produce cutouts in the output folder SNDataResults, and perform forced photometry at the position RAX DECX. However, if the DEC is negative, it should be written as RAX,-DECX; otherwise, argparse reads this as another argument. The log (with all the terminal outputs) will be saved in the data1_path folder SNXXXX (same as the output) as a text file name with SNXXXX.log


---

## 🗂 Data Download

Download Quicklook data (within 7 days):

```bash
python3 subphot_subtract.py -qdl 20250505
```

Download RecentData (within \~30 days):

```bash
python3 subphot_subtract.py -rdl 20250505
```

---

## 🛠 Options and Flags

| Flag       | Description                               |
| ---------- | ----------------------------------------- |
| `-f`       | Folder of images (under `data1_path`)     |
| `-i`       | Individual images to reduce               |
| `-stk`     | Stack images                              |
| `-mp`      | Use multiprocessing (not default)         |
| `-e`       | Email notification (requires SMTP config) |
| `-new`     | Only process new images                   |
| `-b`       | Bands to process (e.g. `-b g r`)          |
| `-sdsscat` | Use SDSS as reference catalog             |
| `-sn`      | Science name(s) (e.g. `-sn SN2023abc`)    |
| `-log`     | Log output to file                        |
| `-o`       | Output folder name                        |
| `-refimg`  | Custom reference image                    |
| `-refcat`  | Custom reference catalog                  |
| `-qdl`     | Download from Quicklook                   |
| `-rdl`     | Download from RecentData                  |
| `-cut`     | Produce image cutouts                     |
| `-un`      | Skip subtraction                          |
| `-up`      | Upload to Fritz                           |
| `-upf`     | Replace existing photometry on Fritz      |
| `-cln`     | Clean intermediate files                  |
| `-mrup`    | Morning round-up subtraction              |
| `-tel`     | Telescope (`LT`, `SEDM`, or `HCT`)        |
| `-lsf`     | List FITS file metadata                   |
| `-fp`      | Forced photometry at specified RA/DEC     |

---

## 📤 Example Output

* Aligned, PSF-matched, and subtracted FITS files.
* Catalogs and photometry tables.
* Optional: uploads to Fritz, image cutouts, and logs.

---

## 📝 Notes

* The system currently uses web scraping to fetch Quicklook/RecentData; more robust methods exist and may be integrated in future updates.
* Some config variables like `sexpath` may be legacy placeholders—keep the path consistent across files to avoid issues.

---

## 📧 Contact

For issues, questions, or contributions, please open an [issue on GitHub](https://github.com/kryanhinds/subphot_pipe/issues) or contact the repository maintainer.

---
