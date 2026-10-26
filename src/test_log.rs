//! Capture of log lines in tests.

use std::io;
use std::sync::{Arc, Mutex};

#[derive(Clone, Default)]
struct Captured(Arc<Mutex<Vec<u8>>>);

impl io::Write for Captured {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Runs `f` and returns its result together with the log lines that it writes at `level` or
/// above.
pub(crate) fn capture<T>(level: tracing::Level, f: impl FnOnce() -> T) -> (T, String) {
    let captured = Captured::default();
    let sink = captured.clone();
    let subscriber = tracing_subscriber::fmt()
        .with_writer(move || sink.clone())
        .with_ansi(false)
        .with_max_level(level)
        .finish();
    let res = tracing::subscriber::with_default(subscriber, f);
    let log = String::from_utf8_lossy(&captured.0.lock().unwrap()).into_owned();
    (res, log)
}
