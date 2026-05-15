# IOBR 2.2.1

## New Features

* **Custom Cache Directory**: Added support for customizing the download cache location via `options(IOBR.cache_dir = "your/path")`. New functions:
  - `get_iobr_cache_dir()` - Get current cache directory path
  - `set_iobr_cache_dir(path)` - Set custom cache directory
  - `reset_iobr_cache_dir()` - Reset to default system cache location
  This enables users to store cached data on shared network drives or any preferred location for easy access across sessions.

## Improvements

- Modified content of some files to adapt to portal display and improve minor errors and prompts (#135).
