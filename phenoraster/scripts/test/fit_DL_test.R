# ----------------------------
# 基于已有拟合参数绘制原始值和平滑DL曲线
debug_pixel_DL_plot <- function(y, t, fit_par, pixel_index, year_sel) {
  cat("原始值 summary:\n")
  print(summary(y))
  
  # 生成密集横轴
  t_dense <- seq(min(t), max(t), by = 1)
  
  # 计算拟合值
  fhat_dense <- doubleLog.Zhang(fit_par, t_dense)
  
  cat("\n拟合值 summary:\n")
  print(summary(fhat_dense))
  
  # 绘图
  plot(t, y, type = "p", col = "blue", pch = 19,
       main = paste("Pixel", pixel_index, "Year", year_sel),
       xlab = "DOY", ylab = "Value")
  lines(t_dense, fhat_dense, col = "red", lwd = 2)
  legend("topleft", legend = c("Original", "DL Fitted"),
         col = c("blue", "red"), pch = c(19, NA), lty = c(NA, 1))
  
  # 返回
  return(list(
    t = t,
    y = y,
    t_dense = t_dense,
    fhat_dense = fhat_dense
  ))
}

# ----------------------------
# 使用示例
# y 和 t 来自你上次成功拟合的像元和年份
# fit_par 来自那次成功拟合得到的参数
res <- debug_pixel_DL_plot(y, t, fit_par, pixel_index = 1222, year_sel = 1995)
