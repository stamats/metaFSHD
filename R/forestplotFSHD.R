## resData: data.frame included in resData_metaFSHD.RData
## ensemblid: single character, ENSEMBL ID
## entrezid: single character, ENTREZ ID
## gene.name: single character, gene name
## print.all: logical, print all matches
## called.in.app: logical, is the function called in an app
forestplotFSHD <- function(resData, ensemblid, entrezid, gene.name, print.all = TRUE,
                           called.in.app = FALSE){
  ## for initialization in app
  if(!missing(gene.name)){
    if(!is.character(gene.name))
      stop("Argument 'gene.name' must be a single *character*!")
    if(length(gene.name) > 1)
      stop("Please provide a *single* gene name!")
    if(!called.in.app & gene.name == "")
      stop("Please provide a *non-empty* gene name!")
    sel <- which(resData$gene_name == gene.name)
    if(length(sel) == 0)
      sel <- grep(gene.name, resData$gene_name)
  }
  if(!missing(entrezid)){
    if(length(entrezid) > 1)
      stop("Please provide a *single* ENTREZ ID!")
    if(!called.in.app & entrezid == "")
      stop("Please provide a *non-empty* ENTREZ ID!")
    sel <- which(resData$entrezid == entrezid)
    if(length(sel) == 0)
      sel <- grep(entrezid, resData$entrezid)
  }
  if(!missing(ensemblid)){
    if(!is.character(ensemblid))
      stop("Argument 'ensemblid' must be a single *character*!")
    if(length(ensemblid) > 1)
      stop("Please provide a *single* ENSEMBL ID!")
    if(!called.in.app & ensemblid == "")
      stop("Please provide a *non-empty* ENSEMBL ID!")
    sel <- which(resData$ENSEMBL == ensemblid)
    if(length(sel) == 0)
      sel <- grep(ensemblid, resData$ENSEMBL)
  }
  if(length(sel) > 1){
    all.sel <- resData[sel,c("gene_name", "ENSEMBL", "entrezid")]
    sel <- sel[1]
    if(!missing(gene.name)){
      if(gene.name != ""){
        warning("Only the first match was used for the provided gene name.")
        print(all.sel)
      } 
    }
    if(!missing(entrezid)){
      if(entrezid != ""){
        warning("Only the first match was used for the provided ENTREZ ID.")
        print(all.sel)
      }
    }
    if(!missing(ensemblid)){
      if(ensemblid != ""){
        warning("Only the first match was used for the provided ENSEMBL ID.")
        print(all.sel)
      }
    }
  }
  if(length(sel) == 0)
    stop("Provided ID not found!")
  
  ind.Mctrl <- grep("Mctrl", names(resData))
  ind.Mfshd <- grep("Mfshd", names(resData))
  ind.SDctrl <- grep("SDctrl", names(resData))
  ind.SDfshd <- grep("SDfshd", names(resData))
  ind.nctrl <- grep("nctrl", names(resData))
  ind.nfshd <- grep("nfshd", names(resData))
  Datasets <- sapply(strsplit(names(resData)[ind.Mfshd], "\\."), "[", 1)
  Datasets[1] <- paste(Datasets[1], "(DUX4)")
  significant <- resData$significant[sel]
  plotData <- data.frame(Dataset = Datasets,
                         mean.fshd = unlist(resData[sel,ind.Mfshd]),
                         mean.ctrl = unlist(resData[sel,ind.Mctrl]),
                         sd.fshd = unlist(resData[sel,ind.SDfshd]),
                         sd.ctrl = unlist(resData[sel,ind.SDctrl]),
                         n.fshd = unlist(resData[sel,ind.nfshd]),
                         n.ctrl = unlist(resData[sel,ind.nctrl]))
  rownames(plotData) <- NULL
  tmp <- escalc(measure = "SMDH", 
                m1i = mean.fshd, m2i = mean.ctrl,
                sd1i = sd.fshd, sd2i = sd.ctrl,
                n1i = n.fshd, n2i = n.ctrl,
                data = plotData, append = TRUE)
  fit <- suppressWarnings(rma(yi, vi, data = tmp, subset = 2:13, measure = "SMDH", method = "REML"))
  FNote <- paste0("(Q = ", formatC(fit$QE, digits=2, format="f"),
                ", df = ", fit$k - fit$p,
                ", p ", metafor:::.pval(fit$QEp, digits=3, showeq=TRUE, sep=" "), "; ",
                "I² = ", formatC(fit$I2, digits=1, format="f"), "%)")
  plotData$mean.diff <- tmp$yi
  z <- qnorm(0.975)
  plotData$ci.low <- tmp$yi - z*sqrt(tmp$vi)
  plotData$ci.upp <- tmp$yi + z*sqrt(tmp$vi)
  plotData$size <- sqrt(tmp$vi)
  plotData$size <- plotData$size/max(plotData$size[-1], na.rm = TRUE)
  plotData$`FSHD` <- paste0(signif(plotData$mean.fshd, 3), " (", 
                            signif(plotData$sd.fshd, 3), ")")
  plotData$`CTRL` <- paste0(signif(plotData$mean.ctrl, 3), " (", 
                            signif(plotData$sd.ctrl, 3), ")")
  sumData <- list(Dataset = "RE Model", 
                  mean.fshd = NA, mean.ctrl = NA, 
                  sd.fshd = NA, sd.ctrl = NA, 
                  n.fshd = NA, n.ctrl = NA, 
                  mean.diff = resData$eff[sel], 
                  ci.low = resData$CI.low[sel], 
                  ci.upp = resData$CI.upp[sel], size = NA, 
                  "FSHD" = "", "CTRL" = "")
  plotData <- rbind(plotData, sumData)
  NA.line <- list(Dataset = "", 
                  mean.fshd = NA, mean.ctrl = NA, 
                  sd.fshd = NA, sd.ctrl = NA, 
                  n.fshd = NA, n.ctrl = NA, 
                  mean.diff = NA, 
                  ci.low = NA, 
                  ci.upp = NA, size = NA, 
                  "FSHD" = "", "CTRL" = "")
  plotData <- rbind(plotData, NA.line)
  plotData$` ` <- paste(rep(" ", 30), collapse = " ")
  plotData$`STD log2-FC (95% CI)` <- paste0(signif(plotData$mean.diff, 3), " (", 
                                           signif(plotData$ci.low, 3), " to ",
                                           signif(plotData$ci.upp, 3), ")")
  plotData$`STD log2-FC (95% CI)`[15] <- ""
  plotData <- plotData[c(1, 15, 2:7, 15, 8:13, 15, 14),]
  plotData$FSHD[plotData$FSHD == "NA (NA)"] <- ""
  plotData$CTRL[plotData$CTRL == "NA (NA)"] <- ""
  plotData$`STD log2-FC (95% CI)`[plotData$`STD log2-FC (95% CI)` == "NA (NA to NA)"] <- ""
  tm <- forest_theme(# Change summary color for filling and borders
                     summary_fill = "black",
                     summary_col = "black",
                     ci_lwd = 1.5,
                     xaxis_lwd = 1.25,
                     refline_lwd = 1.25,
                     arrow_lwd = 1.5,
                     footnote_cex = 0.8)
  TITEL <- "Random-Effects Meta-Analysis of"
  TITEL0 <- paste(strwrap(paste0("ENSEMBL-ID: ", resData$ENSEMBL[sel]), 62), collapse="\n")
  if(resData$gene_name[sel] != ""){
    if(significant){
      TITEL <- paste0(TITEL, " ", resData$gene_name[sel],
                      "*\n", TITEL0)
    }else{
      TITEL <- paste0(TITEL, " ", resData$gene_name[sel],
                      "\n", TITEL0)
    }
  }else{
    TITEL <- paste(TITEL, TITEL0)
  }
  if(resData$entrezid[sel] != ""){
    TITEL1 <- paste(strwrap(paste0("ENTREZ-ID: ", resData$entrezid[sel]), 62), 
                            collapse="\n")
    TITEL <- paste0(TITEL, "\n", TITEL1)
  }
  if(!significant) TITEL <- paste0(TITEL, "\n")
  plotData$size[1] <- 0.5
  plotData$size[nrow(plotData)] <- 0.8
  gg <- forestploter::forest(plotData[,c(1, 12:15)],
                       est = plotData$mean.diff,
                       lower = plotData$ci.low,
                       upper = plotData$ci.upp,
                       sizes = plotData$size,
                       title = TITEL,
                       is_summary = c(rep(FALSE, 16), TRUE),
                       ci_column = 4,
                       arrow_lab = c("FSHD down", "FSHD up"),
                       footnote = "\n\n\n\n\n\nSTD = standardized\nSTD log2-FC = SMD with heteroscedastic variances (SMDH)",
                       theme = tm,
                       xlab = "STD log2-FC")
  gg <- insert_text(gg, text = "Mean (SD) on log2-scale",
                    col = 2:3, part = "header",
                    gp = gpar(fontface = "bold"))
  gg <- add_border(gg, part = "header", row = 1, where = "top")
  gg <- add_border(gg, part = "header", row = 2, where = "bottom")
  gg <- add_border(gg, part = "header", row = 1, col = 2:3, 
                   gp = gpar(lwd = 2))  
  gg <- insert_text(gg,
                    text = "Not included in meta-analysis",
                    row = 2, just = "left",
                    gp = gpar(cex = 0.8))
  gg <- add_border(gg, row = 3, where = "top", gp = gpar(lwd = 1, lty = 2))
  gg <- add_border(gg, row = 17, where = "top", gp = gpar(lwd = 1.5))
  gg <- insert_text(gg, text = FNote,
                    row = 19, just = "left",
                    gp = gpar(cex = 1))
  gg <- insert_text(gg,
                    text = "Microarray",
                    col = 1, row = 4, just = "left",
                    gp = gpar(cex = 1, fontface = "bold"))
  gg <- insert_text(gg,
                    text = "RNA-Seq",
                    col = 1, row = 12, just = "left",
                    gp = gpar(cex = 1, fontface = "bold"))
  if(significant){
    gg <- insert_text(gg, text = "*in 1935 list of Schätzl et al. (2023)",
                      col = 5, row = -1, just = "right",
                      gp = gpar(cex = 0.8, fontface = "bold"))
  }
  row.up <- which(plotData$ci.low > 0)
  if(length(row.up) > 0){
    gg <- edit_plot(gg, col = 4, row = row.up, which = "ci", 
                    gp = gpar(col = "#BC3C29", fill = "#BC3C29", fontface = "plain"))
    row.up[row.up > 1] <- row.up[row.up > 1] + 2
    row.up[row.up > 10] <- row.up[row.up > 10] + 1
    gg <- edit_plot(gg, col = 5, row = row.up, which = "text", 
                    gp = gpar(col = "#BC3C29", fill = "#BC3C29", fontface = "plain"))
  }
  row.dn <- which(plotData$ci.upp < 0)
  if(length(row.dn) > 0){
    gg <- edit_plot(gg, col = 4, row = row.dn, which = "ci", 
                    gp = gpar(col = "#0072B5", fill = "#0072B5", fontface = "plain"))
    row.dn[row.dn > 1] <- row.dn[row.dn > 1] + 2
    row.dn[row.dn > 10] <- row.dn[row.dn > 10] + 1
    gg <- edit_plot(gg, col = 5, row = row.dn, which = "text", 
                    gp = gpar(col = "#0072B5", fill = "#0072B5", fontface = "plain"))
  }
  gg <- edit_plot(gg, row = seq(2,20,by = 2), which = "background", 
            gp = gpar(fill = "#ffffff"))
  gg <- edit_plot(gg, row = seq(1,21,by = 2), which = "background", 
                  gp = gpar(fill = "#EEEEEE"))
  if(length(row.dn) > 0){
    gg <- add_text(gg, text = "*sig. down",
                   col = 5, row = 22, just = "right",
                   gp = gpar(col = "#0072B5", cex = 1, fontface = "italic"),
                   padding = unit(0.5, "cm"))
  }
  if(length(row.up)){
    rw <- ifelse(length(row.dn) > 0, 23, 22)
    gg <- add_text(gg, text = "*sig. up",
                   col = 5, row = rw, just = "right",
                   gp = gpar(col = "#BC3C29", cex = 1, fontface = "italic"),
                   padding = unit(0.5, "cm"))
  }
  gg
}
