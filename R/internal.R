.load_data <- function(name) {
  # 检查缓存
  if (!exists(".iobr_cache", envir = globalenv())) {
    assign(".iobr_cache", new.env(parent = emptyenv()), envir = globalenv())
  }
  cache <- get(".iobr_cache", envir = globalenv())
  
  if (exists(name, envir = cache)) return(get(name, envir = cache))
  
  # 如果是第一次需要 sysdata，加载全部
  if (!exists(".sysdata_loaded", envir = cache)) {
    
    sysdata_path <- system.file("R/sysdata.rda", package = "IOBR", mustWork = FALSE)
    
    # 双重检查：如果没找到，尝试其他方法
    if (!file.exists(sysdata_path)) {
      # 备用方法：直接查找包目录
      pkg_path <- find.package("IOBR", quiet = TRUE)
      if (length(pkg_path) > 0) {
        sysdata_path <- file.path(pkg_path, "R", "sysdata.rda")
      }
    }
    
    if (file.exists(sysdata_path)) {
      load(sysdata_path, envir = cache)
      assign(".sysdata_loaded", TRUE, envir = cache)
      message("Loaded IOBR internal data (", 
              round(file.info(sysdata_path)$size / 1024^2, 1), " MB)")
    } else {
      # 如果还是找不到，设置标记避免重复查找
      assign(".sysdata_loaded", TRUE, envir = cache)
      warning("sysdata.rda not found in IOBR package")
    }
  }
  
  if (exists(name, envir = cache)) return(get(name, envir = cache))
  
  # 尝试 data/ 目录
  tryCatch({
    utils::data(list = name, package = "IOBR", envir = environment())
    data_obj <- get(name, envir = environment())
    assign(name, data_obj, envir = cache)
    return(data_obj)
  }, error = function(e) {
    available <- ls(envir = cache)
    available <- available[!startsWith(available, ".")]
    stop("Data '", name, "' not found in IOBR. Available: ",
         paste(sort(available), collapse = ", "))
  })
}
