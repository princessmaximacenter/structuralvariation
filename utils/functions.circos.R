init_circos_default = function() {
  par(mar = c(1, 1, 1, 1))
  circos.par("start.degree" = 90)
  circos.par("track.height" = 0.15)
  circos.par("cell.padding"=c(0.01,0.01),"track.margin"=c(0.02, 0.02))
  circos.par("canvas.xlim" = c(-1.3, 1.3), "canvas.ylim" = c(-1.3, 1.3))
  circos.initializeWithIdeogram(species = "hg38", plotType = c("ideogram","labels"),labels.cex=1)
}
draw_cna_call = function(plot_seg,draw_loh=F) {
  
  if(draw_loh==T) {
    circos.genomicTrackPlotRegion(plot_seg,ylim = c(0,1),panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = ifelse(value[[1]]=="gain", "red", 
                                                     ifelse(value[[1]] == "loss", "blue",
                                                            ifelse(value[[1]] == "loh","gold","grey"))), 
                         border=NA, ...)
    }) 
  } else {
    plot_seg = plot_seg %>% filter(call!="loh")
    circos.genomicTrackPlotRegion(plot_seg,ylim = c(0,1),panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = ifelse(value[[1]]=="gain", "red", 
                                                     ifelse(value[[1]] == "loss", "blue","grey")), 
                         border=NA, ...)
    })
  }
}
draw_cna_l2fc = function(scaled_seg) {
  
  circos.genomicTrackPlotRegion(scaled_seg,panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = ifelse(value[[1]] > log2foldchange_cutoff, "red", ifelse(value[[1]] < (-log2foldchange_cutoff), "blue","grey")), 
                       border=NA, ytop.column = 1, ybottom = 0, ...)
  })
  
}

draw_cna_maf = function(maf_seg) {
  circos.genomicTrackPlotRegion(maf_seg,panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = "gold",#ifelse(value[[1]] > log2foldchange_cutoff, "red", ifelse(value[[1]] < (-log2foldchange_cutoff), "blue","grey")), 
                       border=NA, ytop.column = 1, ybottom = 0, ...)
  })
}

draw_fusions = function(selected_fusions) {
  if(nrow(selected_fusions)>0){
    sv_start = GRanges(selected_fusions$gup_sv_merged_coordinate) %>% as.data.frame()
    sv_end = GRanges(selected_fusions$gdw_sv_merged_coordinate) %>% as.data.frame()
    
    link_colors = rand_color(nrow(sv_start), transparency = 0.5,luminosity = "bright")
    
    
    circos.genomicLink(sv_start, sv_end, 
                       col=NA,
                       #col = link_colors,
                       border = link_colors)
  }
}

draw_ctx = function(display_ctx) {
  if(nrow(display_ctx)>0){
    sv_start = GRanges(display_ctx$sv_merged_coordinate) %>% as.data.frame()
    sv_end = GRanges(display_ctx$partner_sv_merged_coordinate) %>% as.data.frame()
    
    link_colors = rand_color(nrow(sv_start), transparency = 0.5,luminosity = "bright")
    
    circos.genomicLink(sv_start, sv_end, 
                       col=NA,
                       border = link_colors)
  }
}
draw_intra_svs = function(display_intra_svs) {
  if(nrow(display_intra_svs)>0){
    intra_svs_start = display_intra_svs[,c("seqnames","start")]
    intra_svs_start$end = intra_svs_start$start+1
    intra_svs_end = display_intra_svs[,c("seqnames","end")]
    intra_svs_end$start = intra_svs_start$end-1
    
    #link_colors = rand_color(nrow(intra_svs_start), transparency = 0.5,luminosity = "bright")
    
    circos.genomicLink(intra_svs_start, intra_svs_end,
                       col=add_transparency(display_intra_svs$color, 0.5),
                       border = add_transparency(display_intra_svs$color, 0.5))
  }
}
