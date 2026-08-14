use thiserror::Error;

#[derive(Error, Debug)]
pub enum Error {
    #[error("Fitsio error: {0}")]
    FitsioError(#[from] fitsio::errors::Error),
    #[error("file not found: {0}")]
    FileNotFound(String),
    #[error("invalid data: {0}")]
    InvalidData(String),
    #[error("time reference invalid: {0}")]
    TimeReferenceInvalid(String),
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("unknown detector: {0}")]
    UnknownDetector(String),
    #[error("unknown error occurred")]
    Unknown,
}
