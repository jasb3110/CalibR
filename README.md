# CalibR

## Contents
- [Proposal](#proposal)
- [Instruction](#instruction)
- [R code](#r-code)
- [Reference](#reference)

## Proposal 

"calibR4" is an R function designed to replicate and enhance radiocarbon (C-14) calibration, inspired by Calib 8.02, providing:

- Flexible calibration curves (intcal20, marine20, shcal20, etc.), including mixed marine–terrestrial cases.
- Explicit sigma control (1σ/2σ/3σ or any probability between 0–1).Reproducible outputs: native calibration PNG, ggplot2-enhanced PNG, CSV with probability density, and input/output tables as PNG.
- Smart labels for mean, median, maximum values with ggrepel.

Below, I’ve included the script and an example using published data from Guiñez et al., 2014[(3)](#reference).
The input table is displayed below. It is crucial to sort the data in this specific format.

|[![Table 1.](calout/mejillones.input.png)](https://github.com/jasb3110/CalibR/blob/45c5c0bdf3854b9ce31530ca285ea5a6f03fc1c9/calout/mejillones.input.png)|
|:--:| 
|*Table. Input table of Radiocarbon samples pulled out Mejillones core (Guiñez et al. 2014)*|

In the first plot, you can refer to the format of the outcome. The x-axis represents the range of calibrated ages, while the y-axis shows the density probability of calibration. The black, green, and gray lines indicate the maximum probability, median, and mean calibrated ages, respectively. Additionally, the blue and red lines denote the calibration limits according to one standard deviation (one sigma: 68%). You can adjust the sigma value if desired.

|[![Plot.](calout/mejillones/mejillones-mejillones-18%20sample%20at%2054%20cm.png)](https://github.com/jasb3110/CalibR/blob/45c5c0bdf3854b9ce31530ca285ea5a6f03fc1c9/calout/mejillones/mejillones-mejillones-18%20sample%20at%2054%20cm.png)|
|:--:| 
|*Picture 1. Outcome plot of sample on 54cm colour-scale (Guiñez et al. 2014)*|

The second plot is similar to the previous one; however, it is presented in grayscale.

|[![Plot.](calout/mejillones-gray%20version/mejillones-mejillones-18%20sample%20at%2054%20cm.png)](https://github.com/jasb3110/CalibR/blob/45dfc25c6e44d5bb47df28493aa66481c0ccbc4a/calout/mejillones-gray%20version/mejillones-mejillones-18%20sample%20at%2054%20cm.png)|
|:--:| 
|*Picture 2. Outcome plot of sample on 54cm in gray-scale (Guiñez et al. 2014)*|

Finally, This script is created an outcome table where you are able to find the maximum, median, and mean calibrate age in columns. the outcome table will be saved together with plots [(Table 2.)](calout/mejillones/mejillones.output.png) and I attached a function in source [(4)](https://github.com/jasb3110/CalibR/blob/36366a16fc9cf5e4c5f070f9b17be2f357915dc5/calib.R).

|[![Table 2.](calout/mejillones.output.png)](https://github.com/jasb3110/CalibR/blob/45dfc25c6e44d5bb47df28493aa66481c0ccbc4a/calout/mejillones.output.png)|
|:--:| 
|*Table 2. Output table of Radiocarbon samples pulled out Mejillones core (Guiñez et al. 2014)*|

## Instruction

1. Select the whole script and pulse Ctrl+ Enter.
2. To wait for its outcomes when you will hear Mario Bross sound the means it is finished. 
3. Bon appetit!!

### WARNING:

Some plots may appear distorted, so you should adjust this script to improve the margins, resolution, or font of the images as needed.

## R code

Finally, it were showed a source of Calib4.R[(4)](#reference). 
```markdown
#########################################################################
calibR4<- function(
    input=input,
    sigma=c(1,2,3,.68,.95,.99,"1s","2s","3s","1sigma","2sigma","3sigma"),
    curve=c(1,2,3,"intcal20","marine20","shcal20","marine13","shcal13",
            "nh1","nh2","nh3","sh1-2","sh3","nh1_monthly","nh1_monthly","nh2_monthly",
            "nh3_monthly","sh1-2_monthly","sh3_monthly","kure","levinKromer","santos"),
    colour=c("default","color","colour",1,TRUE,"yes","minimal","gray","grey",0,FALSE,"no"),
    show.table=c(TRUE,1,"yes","default",FALSE,0,"no"),
    show.plot=c("both",TRUE,1,"yes","minimal","default",FALSE,0,"no"),
    width.native = 200, height.native  = 200, units.native  = 'mm', res.native = 1200,
    dpi.calib = 900, width.calib = 250, height.calib = 159, units.calib="mm",
    round = 2,
    csv_semicolon = FALSE  # TRUE para CSV con ';'
){
  ##############################################################################
  # License: GNU
  # Autor: José Solís, november 2025
  # email: solisbenites.jose@gmail.com
  ##############################################################################
  # Paquetes
  req <- c("IntCal","ggplot2","ggh4x","gridExtra","magrittr","scales","ggrepel","dplyr","progress")
  ok <- sapply(req, require, character.only=TRUE, quietly=TRUE)
  if (!all(ok)) stop("Missing packages: ", paste(req[!ok], collapse=", "))
  has_beep <- requireNamespace("beepr", quietly=TRUE)
  beep <- function(i){ if (has_beep) beepr::beep(i) }
  
  begin <- Sys.time()
  if (missing(input)) { beep(7); stop("Missing 'input'") }
  
  # ------------------------------ Normalizadores ------------------------------
  to_bool <- function(x, yes=c("both","yes","true","1","t","y"), no=c("no","false","0","f","n")){
    if (is.logical(x)) return(x)
    if (is.numeric(x)) return(x != 0)
    xx <- tolower(as.character(x))
    if (xx %in% yes)  return(TRUE)
    if (xx %in% no)   return(FALSE)
    NA
  }
  norm_colour <- function(x){
    xx <- tolower(as.character(x))
    if (xx %in% c("default","color","colour","1","yes","true")) return("color")
    if (xx %in% c("minimal","gray","grey","0","no","false"))    return("gray")
    "gray"
  }
  norm_showplot <- function(x){
    xx <- tolower(as.character(x))
    if (xx %in% c("both","yes","true","1")) return("both")
    if (xx %in% c("minimal","default"))      return("file")  # solo guardar a disco
    if (xx %in% c("no","false","0"))         return("none")
    "file"
  }
  map_sigma <- function(s){
    if (length(s) != 1) stop("Sigma must be a single value.")
    if (is.numeric(s)){
      if (dplyr::between(s,0,1)) return(s)
      if (s %in% c(1,2,3)) return(c(.68,.95,.99)[match(s, c(1,2,3))])
      stop("Invalid numeric sigma: use 0-1 (prob) or 1/2/3.")
    } else {
      ss <- tolower(as.character(s))
      look  <- c("1s","2s","3s","1sigma","2sigma","3sigma")
      vals  <- c(.68,.95,.99,.68,.95,.99)
      if (ss %in% look) return(vals[match(ss,look)])
      if (grepl("^0?\\.\\d+$", ss)) return(as.numeric(ss))
      stop("Invalid sigma string. Use '1s','2s','3s','1sigma', etc.")
    }
  }
  
  # ------------------------------- Datos de entrada ---------------------------
  dd  <- as.data.frame(input[-1,], stringsAsFactors = FALSE)
  dct <- getwd()
  
  # Carpeta raíz de salida
  out_root <- file.path(dct, "calout")
  if (!dir.exists(out_root)) { dir.create(out_root, recursive=TRUE); message("To create 'calout' file") }
  
  sanitize_name <- function(x){
    x <- as.character(x); x[is.na(x)] <- ""
    x <- gsub("[^A-Za-z0-9._-]+","_",x); x[nchar(x)==0] <- "LAB_UNNAMED"; x
  }
  get_lab_dir <- function(lab_value){
    lab_dir <- file.path(out_root, sanitize_name(lab_value))
    if (!dir.exists(lab_dir)) { dir.create(lab_dir, recursive=TRUE); message("Created folder: ", lab_dir) }
    lab_dir
  }
  
  # Nombres esperados
  namess <- c("Lab","Sample","X14C.BP","X14C.Age.SD","Lab.Error","Age.Span",
              "Uncorrected.14C","Uncorrected.14C.SD","d13C","d13C.SD",
              "Delta.R","Delta.R.SD","Marine.Carbon","Description","CalCurve","Depth")
  if (ncol(dd) != length(namess) || !identical(colnames(dd), namess)){
    beep(7); stop("Input columns do not exactly match the expected schema.")
  }
  
  # Flags y sigma
  colour_mode <- norm_colour(colour[1])
  plot_mode   <- norm_showplot(show.plot[1])
  table_plot  <- to_bool(show.table[1]); if (is.na(table_plot)) { beep(7); stop("Invalid 'show.table'.") }
  sigma_prob  <- map_sigma(sigma[1])
  use_color   <- identical(colour_mode, "color")
  
  # Campos de salida
  add_cols <- c("mean","lower","upper","median","max","error","percent")
  for (v in add_cols) dd[[v]] <- NA_real_
  
  # ----------------------------- Progreso -------------------------------------
  pb <- progress_bar$new(
    format = "(:spin) [:bar] :percent [Elapsed: :elapsedfull | ETA: :eta]",
    total = nrow(dd),
    complete = "=",
    incomplete = "-",
    current = ">",
    clear = FALSE,
    width = 100
  )
  
  # ------------------------------ Calibración ---------------------------------
  beep(2); beep(2); beep(2)
  
  for (i in seq_len(nrow(dd))) {
    pb$tick()
    
    label.name <- if (isTRUE(is.na(dd$Depth[i]) || dd$Depth[i]==""))
      paste0(dd$Lab[i],"-",dd$Sample[i],"-",i) else
        paste0(dd$Lab[i],"-",dd$Sample[i]," sample at ",dd$Depth[i]," cm")
    
    rsv   <- suppressWarnings(as.numeric(dd$Delta.R[i]))
    sdrsv <- suppressWarnings(as.numeric(dd$Delta.R.SD[i]))
    c14   <- suppressWarnings(as.numeric(dd$X14C.BP[i]))
    sdc14 <- suppressWarnings(as.numeric(dd$X14C.Age.SD[i]))
    marC  <- suppressWarnings(as.numeric(dd$Marine.Carbon[i]))
    
    if (is.na(dd$CalCurve[i]) || length(dd$CalCurve[i])!=1){
      beep(7); stop("You must choose a calibration curve for: ", label.name)
    }
    
    # Curva
    curve_i <- dd$CalCurve[i]
    if (is.numeric(curve_i)) {
      curve_i <- c("intcal20","marine20","shcal20")[as.integer(curve_i)]
    }
    curve_i <- tolower(as.character(curve_i))
    
    # Mezcla marino/terrestre
    if (is.na(marC) || !dplyr::between(marC,0,100)){
      beep(7); stop("Marine.Carbon must be between 0 and 100 (", label.name,")")
    }
    if (marC == 100) curve_i <- "marine20"
    if (marC == 0){
      if (is.na(rsv))   rsv <- 0
      if (is.na(sdrsv)) sdrsv <- 0
    }
    if (marC > 0 && marC < 100){
      cc2 <- if (curve_i %in% c("intcal20","marine20","intcal13","nh1","nh2","nh3",
                                "nh1_monthly","nh2_monthly","nh3_monthly")) "IntCal20" else "SHCal20"
      try(IntCal::mix.ccurves(proportion = marC/100, cc1="Marine20", cc2=cc2,
                              name = "mixed.marine.terrigenous.14C",
                              dir = c(), offset = c(0,0), sep = "\t"),
          silent=TRUE)
      curve_i <- "mixed.marine.terrigenous"
    }
    
    # Guardas NA
    if (is.na(c14) || is.na(sdc14)){
      warning("Cannot calibrate (NA) at ", label.name); next()
    }
    if (is.na(rsv))   { rsv <- 0;   warning("\u0394R is NA; using 0 in ", label.name) }
    if (is.na(sdrsv)) { sdrsv <- 0; warning("\u0394R.SD is NA; using 0 in ", label.name) }
    
    # Rango de curva
    curv <- IntCal::ccurve(curve_i)
    pp <- sum((c14 - sdc14 - rsv - sdrsv) <= min(curv$V2 - curv$V3),
              (c14 + sdc14 + rsv + sdrsv) >= max(curv$V2), na.rm=TRUE)
    
    # Carpeta por Lab
    lab_dir <- get_lab_dir(dd$Lab[i])
    
    # Fuera de rango => PNG nativo vacío y continuar
    if (pp == 1){
      plot.native <-
        ggplot2::ggplot()+
        ggplot2::theme_void()+
        ggplot2::geom_text(ggplot2::aes(0,0,label=paste0("Cannot calibrate: ",label.name)))
      
      ggplot2::ggsave(filename=file.path(lab_dir, paste0(label.name,".native.png")),
                      dpi=res.native, width=width.native, height=height.native,
                      units=units.native, plot=plot.native)
      message("\nOut of curve range: ", label.name)
      next()
    }
    
    # Calibración
    res <- IntCal::calibrate(age=c14, error=sdc14, cc=curve_i,
                             prob=sigma_prob, yr.steps=1, threshold=5e-04,
                             rounded=round, reservoir=c(rsv,sdrsv), legend.cex=1)
    
    # Guardar PNG nativo (base plot del paquete)
    png(file.path(lab_dir, paste0(label.name,".native.png")),
        width=width.native, height=height.native, units=units.native, res=res.native)
    invisible(IntCal::calibrate(age=c14, error=sdc14, cc=curve_i,
                                prob=sigma_prob, yr.steps=1, threshold=5e-04,
                                rounded=round, reservoir=c(rsv,sdrsv), legend.cex=1))
    dev.off()
    
    if (plot_mode == "both" && interactive()){
      print(res)
    }
    
    # Extraer densidad e intervalos
    gg <- as.data.frame(res[[2]])
    if (isTRUE(csv_semicolon)) {
      utils::write.csv2(gg, file=file.path(lab_dir, paste0(label.name,".",sigma_prob,".probability.csv")),
                        row.names=FALSE)
    } else {
      utils::write.csv(gg,  file=file.path(lab_dir, paste0(label.name,".",sigma_prob,".probability.csv")),
                       row.names=FALSE)
    }
    
    dr  <- as.data.frame(res[[1]]); colnames(dr) <- c("cal.BP","V2")
    den <- dr
    
    # Intervalo principal
    if (nrow(gg) == 1){
      dd$lower[i]   <- gg$from[1]; dd$upper[i] <- gg$to[1]; dd$percent[i] <- gg$perc[1]
    } else {
      j <- which.max(gg$perc)
      dd$lower[i]   <- gg$from[j]; dd$upper[i] <- gg$to[j]; dd$percent[i] <- gg$perc[j]
    }
    
    # Mediana, media y modo (máximo)
    sel <- which(den$cal.BP >= dd$upper[i] & den$cal.BP <= dd$lower[i])
    if (length(sel) > 1){
      dens_sel <- den$V2[sel]; x_sel <- den$cal.BP[sel]
      csum <- cumsum(dens_sel)/sum(dens_sel)
      dd$median[i] <- x_sel[which.max(csum >= 0.5)]
      dd$mean[i]   <- sum(x_sel * dens_sel) / sum(dens_sel)
      dd$max[i]    <- x_sel[which.max(dens_sel)]
      dd$error[i]  <- max(abs(dd$lower[i]-dd$max[i]), abs(dd$max[i]-dd$upper[i]))
    }
    
    # Asegurar límites 0/rango
    limx <- range(IntCal::ccurve(curve_i)$V1, na.rm=TRUE)
    dd$max[i]   <- min(max(dd$max[i], 0), limx[2]*.9999)
    dd$lower[i] <- min(dd$lower[i],      limx[2]*.9999)
    dd$upper[i] <- max(dd$upper[i],      0)
    if (!is.na(dd$mean[i])){
      dd$mean[i] <- min(max(dd$mean[i],0), limx[2]*.9999)
    }
    
    # Polígonos y tabla ordenada/sin duplicados
    dr2 <- den[is.finite(den$cal.BP) & is.finite(den$V2) & den$V2 >= 0, ]
    dr2 <- dr2[order(dr2$cal.BP), ]
    dr2 <- dplyr::distinct(dr2, cal.BP, .keep_all=TRUE)
    x_left  <- min(dr2$cal.BP); x_right <- max(dr2$cal.BP)
    poly_dr <- rbind(
      c(x_left,  0),
      c(x_left,  dr2$V2[1]),
      as.matrix(dr2[,c("cal.BP","V2")]),
      c(x_right, 0),
      c(x_left,  0)
    )
    poly_dr <- as.data.frame(poly_dr); names(poly_dr) <- c("x","y")
    
    # Área de confianza
    da <- dr2[dr2$cal.BP > dd$upper[i] & dr2$cal.BP < dd$lower[i], c("cal.BP","V2")]
    if (nrow(da) == 0) warning("Empty confidence polygon in ", label.name)
    names(da) <- c("x","y")
    
    # Límites y etiquetas
    xr       <- range(poly_dr$x, na.rm=TRUE); pad <- diff(xr)*0.05
    xlim_rev <- c(xr[2]+pad, xr[1]-pad)
    ylim_ok  <- c(0, max(poly_dr$y, na.rm=TRUE)*1.02)
    
    qx   <- stats::quantile(dr2$cal.BP, na.rm=TRUE)
    qyM  <- max(dr2$V2, na.rm=TRUE)
    lab_x <- mean(c(qx[4], qx[5]))
    lab_y <- 0.96*qyM
    
    # Altura y para x dada
    y_at <- function(x) {
      if (is.na(x)) return(NA_real_)
      if (abs(x) < 1e-12) {
        y0 <- suppressWarnings(max(den$V2[den$cal.BP == 0], na.rm = TRUE))
        if (is.finite(y0)) return(y0)
      }
      approx(dr2$cal.BP, dr2$V2, x, rule = 2)$y
    }
    
    cols <- if (use_color) {
      c(Mean = "gray50", Median = "green",  Maximum = "black")
    } else {
      c(Mean = "white", Median = "white", Maximum = "white")
    }
    
    eti <- data.frame(
      x   = c(dd$mean[i],   dd$median[i],   dd$max[i]),
      y   = c(y_at(dd$mean[i]), y_at(dd$median[i]), y_at(dd$max[i])),
      lab = factor(c("Mean","Median","Maximum"), levels = names(cols)),
      stringsAsFactors = FALSE
    )
    eti <- stats::na.omit(eti)
    
    # ------------------------------ Gráfico final ------------------------------
    plotting <-
      ggplot2::ggplot(dr2, ggplot2::aes(x = cal.BP, y = V2)) +
      ggplot2::geom_polygon(
        data = poly_dr,
        ggplot2::aes(x = x, y = y, group = 1),
        fill = "white", colour = "black", linewidth = .1, inherit.aes = FALSE
      ) +
      ggplot2::geom_area(
        data = da, ggplot2::aes(x = x, y = y),
        fill = "gray90", colour = "black", linewidth = .1, inherit.aes = FALSE
      ) +
      ggplot2::coord_cartesian(ylim = ylim_ok, expand = FALSE) +
      ggplot2::geom_segment(
        ggplot2::aes(x = dd$mean[i], xend = dd$mean[i], y = 0, yend = y_at(dd$mean[i])),
        color = if (use_color) "gray50" else "gray60", alpha = .6, linewidth = .2
      ) +
      ggplot2::geom_segment(
        ggplot2::aes(x = dd$median[i], xend = dd$median[i], y = 0, yend = y_at(dd$median[i])),
        color = if (use_color) "green" else "gray25", alpha = 1, linewidth = .35
      ) +
      ggplot2::geom_segment(
        ggplot2::aes(x = dd$max[i], xend = dd$max[i], y = 0, yend = y_at(dd$max[i])),
        color = "black", alpha = 1, linewidth = .5
      ) +
      ggplot2::scale_x_reverse(
        limits = xlim_rev,
        guide  = ggplot2::guide_axis(minor.ticks = TRUE),
        breaks = scales::pretty_breaks(n = 4),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::scale_y_continuous(
        guide  = ggplot2::guide_axis(minor.ticks = TRUE),
        breaks = scales::pretty_breaks(n = 4),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::annotate("text", x = lab_x, y = lab_y,
                        label = paste0("Cal. age: ", dd$max[i], " \u00B1 ", dd$error[i]),
                        size = 4) +
      ggplot2::annotate("text", x = lab_x, y = lab_y*0.97,
                        label = paste0("\u0394R = ", rsv, " \u00B1 ", sdrsv), size = 4) +
      ggplot2::annotate("text", x = lab_x, y = lab_y*0.94,
                        label = paste0("Probability: ", base::round(dd$percent[i], 2), "%"), size = 4) +
      ggplot2::annotate("text", x = lab_x, y = lab_y*0.91,
                        label = paste0("Cal. curve: ", curve_i), size = 4) +
      ggrepel::geom_text_repel(
        data= eti, 
        ggplot2::aes(x=x,y=y,label=lab,colour=lab,segment.colour="black",alpha = 0.5),
        bg.color = if (use_color) scales::alpha("white", 0.5) else scales::alpha("grey30", 0.5),
        bg.r = 0.05, force = 1, vjust = 0, hjust = 0,
        min.segment.length = grid::unit(.15,'lines'),
        nudge_y = 0, nudge_x = 0,
        fontface = 'bold', label.r = .3, point.padding = grid::unit(.5,'lines'),
        segment.ncp = 3, segment.size  = 0.05,
        box.padding = grid::unit(1.25,'lines'),
        arrow = grid::arrow(ends = "last", type = "closed",angle = 15, length = grid::unit(.1, "lines")),
        max.iter = Inf, show.legend = FALSE, verbose=FALSE, inherit.aes = TRUE
      )+
      ggplot2::scale_colour_manual(values = cols, guide = "none") +
      ggplot2::labs(
        title = "Relative probability of sample",
        x = "Cal yr BP",
        y = expression(paste("Density"))
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.ticks.length = grid::unit(0.2, "cm"),
        ggh4x.axis.ticks.length.minor = ggplot2::rel(0.5),
        axis.ticks  = ggplot2::element_line(linewidth = 0.5),
        axis.text.x = ggplot2::element_text(size = 11, colour = "black", face = "bold", vjust = 0),
        axis.text.y = ggplot2::element_text(size = 11, colour = "black", face = "bold", hjust = 1),
        axis.title  = ggplot2::element_text(size = 14, face = "bold"),
        title       = ggplot2::element_text(size = 16, colour = "black", face = "bold")
      )
    
    # Guardar figura final
    ggplot2::ggsave(filename=file.path(lab_dir, paste0(label.name,".png")),
                    dpi=dpi.calib, width=width.calib, height=height.calib,
                    units=units.calib, plot=plotting)
    
    if (plot_mode == "both" && interactive()){
      print(plotting)
    }
  } # fin for
  
  # ------------------------------ Tablas de salida -----------------------------
  dd$mean    <- base::round(dd$mean)
  dd$percent <- base::round(dd$percent)
  dd$max[which(dd$max==1)] <- 0
  
  # Etiquetas bonitas para tabla de entrada
  l <- c("Lab code", "Sample code",expression(phantom()^14*C~"\n(yrs BP)"),
         expression(phantom()^14*C~SD~"\n(yrs BP)"), "Lab Error", "Age Span",
         "Uncorrected~delta^14*C","Uncorrected~SD~delta^14*C","delta^13*C",
         "delta^13*C~SD" ,"Delta*R~(yrs)",
         "Delta*R~SD~(yrs)","Marine carbon\n(%)","Description","Calibration~curve",
         "Depth~(cm)")
  d1 <- input[-1,]; colnames(d1) <- l
  
  myt <- gridExtra::ttheme_minimal(
    base_size=12, base_colour="black", parse=TRUE,
    padding=grid::unit(c(2,2), "mm"),
    colhead=list(fg_params=list(parse=TRUE), fontface=4L, bg_params=list(fill="gray90"))
  )
  
  # --- Guardado por número de Labs ---
  labs <- unique(dd$Lab)
  
  if (length(labs) == 1) {
    # Un solo Lab: todo junto en su carpeta
    lab <- labs[1]
    lab_dir <- get_lab_dir(lab)
    
    png(file.path(lab_dir, paste0(sanitize_name(lab),".input.png")),
        width = 20+ncol(d1)*425/15, height = 20+100/19*nrow(d1), units='mm', res=1200)
    gridExtra::grid.table(d1, rows=NULL, theme=myt); dev.off()
    
    png(file.path(lab_dir, paste0(sanitize_name(lab),".output.png")),
        width = 20+ncol(dd)*525/22, height = 20+100/19*nrow(dd), units='mm', res=1200)
    gridExtra::grid.table(dd, rows=NULL, theme=myt); dev.off()
    
    if (isTRUE(csv_semicolon)){
      utils::write.csv2(dd, file=file.path(lab_dir, paste0(sanitize_name(lab),".calibrated.csv")), row.names=FALSE)
    } else {
      utils::write.csv(dd,  file=file.path(lab_dir, paste0(sanitize_name(lab),".calibrated.csv")),  row.names=FALSE)
    }
    
    if (isTRUE(table_plot) && interactive()){
      gridExtra::grid.table(d1, rows=NULL, theme=myt)
      gridExtra::grid.table(dd, rows=NULL, theme=myt)
    }
    
  } else {
    # ≥2 Labs: archivos por Lab + CONSOLIDADO GLOBAL con nombres concatenados
    for (lab in labs) {
      lab_dir <- get_lab_dir(lab)
      
      d1_lab <- d1[ d1$`Lab code` == lab, , drop=FALSE ]
      dd_lab <- dd[ dd$Lab       == lab, , drop=FALSE ]
      
      if (nrow(d1_lab) > 0) {
        png(file.path(lab_dir, paste0(sanitize_name(lab),".input.png")),
            width = 20+ncol(d1_lab)*425/15, height = 20+100/19*nrow(d1_lab), units='mm', res=1200)
        gridExtra::grid.table(d1_lab, rows=NULL, theme=myt); dev.off()
      }
      if (nrow(dd_lab) > 0) {
        png(file.path(lab_dir, paste0(sanitize_name(lab),".output.png")),
            width = 20+ncol(dd_lab)*525/22, height = 20+100/19*nrow(dd_lab), units='mm', res=1200)
        gridExtra::grid.table(dd_lab, rows=NULL, theme=myt); dev.off()
      }
      
      if (nrow(d1_lab) > 0) {
        if (isTRUE(csv_semicolon)){
          utils::write.csv2(d1_lab, file=file.path(lab_dir, paste0(sanitize_name(lab),".calibrated.csv")), row.names=FALSE)
        } else {
          utils::write.csv(d1_lab,  file=file.path(lab_dir, paste0(sanitize_name(lab),".calibrated.csv")),  row.names=FALSE)
        }
      }
    }
    
    # -------- CONSOLIDADO GLOBAL con prefijo = concat(Lab codes sanitizados) --------
    cons_prefix <- paste0(sanitize_name(labs), collapse = "_")
    
    png(file.path(out_root, paste0(cons_prefix, ".input.png")),
        width = 20+ncol(d1)*425/15, height = 20+100/19*nrow(d1), units='mm', res=1200)
    gridExtra::grid.table(d1, rows=NULL, theme=myt); dev.off()
    
    png(file.path(out_root, paste0(cons_prefix, ".output.png")),
        width = 20+ncol(dd)*525/22, height = 20+100/19*nrow(dd), units='mm', res=1200)
    gridExtra::grid.table(dd, rows=NULL, theme=myt); dev.off()
    
    if (isTRUE(csv_semicolon)){
      utils::write.csv2(d1, file=file.path(out_root, paste0(cons_prefix, ".input.csv")),  row.names=FALSE)
      utils::write.csv2(dd, file=file.path(out_root, paste0(cons_prefix, ".output.csv")), row.names=FALSE)
    } else {
      utils::write.csv(d1, file=file.path(out_root, paste0(cons_prefix, ".input.csv")),  row.names=FALSE)
      utils::write.csv(dd, file=file.path(out_root, paste0(cons_prefix, ".output.csv")), row.names=FALSE)
    }
  }
  
  end <- Sys.time()
  message("Working time: ", base::round(as.numeric(difftime(end, begin, units="mins")), 2)," minutes")
  beep(8)
  invisible(d1)
}
################################################################################
```
## Reference

  - Intcal package in R (2022). https://cran.r-project.org/web/packages/IntCal/index.html
  - Stuiver, M., & Reimer, P. J. (1993). EXTENDED 14C DATA BASE AND REVISED CALIB 3.014C AGE CALIBRATION PROGRAM. Radiocarbon, 35(1), 215–230. https://doi.org/10.14210/bjast.v17.n2.pNB5-8
  - Guiñez, M., Valdés, J., Sifeddine, A., Boussafir, M., & Dávila, P. M. (2014). Anchovy population and ocean-climatic fluctuations in the Humboldt Current System during the last 700 years and their implications. Palaeogeography, Palaeoclimatology, Palaeoecology, 415, 210–224. https://doi.org/10.1016/j.palaeo.2014.08.026
  - https://github.com/jasb3110/CalibR/blob/02cb2fef4994c204081e9a7e28c4e1e55471d8ce/calib.R  - version 1
  - https://github.com/jasb3110/CalibR/blob/4852cd55965a5c007e0de92ed9fe4c5273570682/calib2.R - version 2
  - https://github.com/jasb3110/CalibR/blob/45dfc25c6e44d5bb47df28493aa66481c0ccbc4a/calib3.R - version 3
  - https://github.com/jasb3110/CalibR/blob/e8c32d088cd90f3215f4615d3bac1863759617ab/calib4.R - version 4 - Last version
