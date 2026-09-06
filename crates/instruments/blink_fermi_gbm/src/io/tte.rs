/*
glg_tte_n0_190101_00z_v00.fit.gz
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU      30   ()
  1  EBOUNDS       1 BinTableHDU     51   128R x 3C   [1I, 1E, 1E]
  2  EVENTS        1 BinTableHDU     38   4007335R x 2C   [1D, 1I]
  3  GTI           1 BinTableHDU     39   1R x 2C   [1D, 1D]
*/

use crate::types::Detector;

/// 一个探头一小时的 TTE。文件是 gzip 压缩的，cfitsio 透明解压。
pub struct TteFile {
    pub detector: Detector,
    time: Vec<f64>,
    pha: Vec<i16>,
    gti_start: Vec<f64>,
    gti_stop: Vec<f64>,
}

impl TteFile {
    pub fn from_fits_file(path: &str, detector: Detector) -> Result<Self, fitsio::errors::Error> {
        let mut fptr = fitsio::FitsFile::open(path)?;

        let events = fptr.hdu("EVENTS")?;
        let time = events.read_col::<f64>(&mut fptr, "TIME")?;
        let pha = events.read_col::<i16>(&mut fptr, "PHA")?;

        let gti = fptr.hdu("GTI")?;
        let gti_start = gti.read_col::<f64>(&mut fptr, "START")?;
        let gti_stop = gti.read_col::<f64>(&mut fptr, "STOP")?;

        Ok(Self {
            detector,
            time,
            pha,
            gti_start,
            gti_stop,
        })
    }

    pub fn len(&self) -> usize {
        self.time.len()
    }

    pub fn is_empty(&self) -> bool {
        self.time.is_empty()
    }

    pub fn time(&self) -> &[f64] {
        &self.time
    }

    pub fn pha(&self) -> &[i16] {
        &self.pha
    }

    /// GTI 各行 (START, STOP)。逐小时文件通常只有一行，在 SAA 进入处截止。
    pub fn gti_rows(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        self.gti_start
            .iter()
            .zip(self.gti_stop.iter())
            .map(|(a, b)| (*a, *b))
    }

    /// 本文件 GTI 与 `[from, to]` 的交集长度（秒）。
    pub fn gti_seconds_within(&self, from: f64, to: f64) -> f64 {
        self.gti_start
            .iter()
            .zip(self.gti_stop.iter())
            .map(|(start, stop)| (stop.min(to) - start.max(from)).max(0.0))
            .sum()
    }

    /// 时间回跳处数。事例流有序是 k 路归并与搜索窗长判据的前提；
    /// SVOM/GRM 上出现过整段数据写两遍导致的回跳，这里一并记账。
    pub fn time_reversals(&self) -> usize {
        self.time.windows(2).filter(|pair| pair[1] < pair[0]).count()
    }
}
