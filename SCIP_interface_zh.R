library(shiny)
library(DT)
library(plotrix)
options(warn=-1)
options(encoding = "UTF-8")

Sys.setlocale("LC_ALL", "en_US.UTF-8")

# https://stackoverflow.com/questions/43635846/
BinMean <- function (vec, every, na.rm = TRUE) {
  n <- length(vec)
  x <- .colMeans(vec, every, n %/% every, na.rm)
  r <- n %% every
  if (r) x <- c(x, mean.default(vec[(n - r + 1):n], na.rm = na.rm))
  x
}

# config_file=read.table("./interface_config.txt",comment.char="",sep="\t",quote="")
# temp_file_dir=config_file[2,2]
# username=config_file[4,2]
# regions_list=read.table(paste(temp_file_dir,"/../user_data/",config_file[1,2],".pipeline_summary.txt",sep=""),comment.char="",sep="\t",quote="")

# https://stackoverflow.com/questions/63882483
mymodal <- function (..., title = NULL, footer = modalButton("Dismiss"), 
                     size = c("m", "s", "l"), easyClose = FALSE, fade = TRUE, idcss = "") 
{
  size <- match.arg(size)
  cls <- if (fade) 
    "modal fade"
  else "modal"
  div(id = "shiny-modal", class = cls, tabindex = "-1", `data-backdrop` = if (!easyClose) 
    "static", `data-keyboard` = if (!easyClose) 
      "false", div(class = paste("modal-dialog", idcss), class = switch(size, 
                                                                        s = "modal-sm", 
                                                                        m = NULL, 
                                                                        l = "modal-lg"), 
                   div(class = "modal-content", 
                       if (!is.null(title)) 
                         div(class = "modal-header", tags$h4(class = "modal-title", 
                                                             title)
                         ), 
                       div(class = "modal-body", ...), 
                       if (!is.null(footer)) 
                         div(class = "modal-footer", footer))
      ), 
    tags$script("$('#shiny-modal').modal().focus();"))
}

# sort by priority, then by size
# regions_name_sort1=regions_list[sort(regions_list$V5-regions_list$V4+1,index.return=T,decreasing=T)$ix,]
# regions_name=regions_name_sort1$V17[sort(regions_name_sort1$V18,index.return=T)$ix]
# results_out=paste(config_file[3,2],"/SCIP_results.txt",sep="")
# if (file.exists(results_out)==F){
#   file.create(results_out)
# }

ui=fluidPage(
  tags$head(
    tags$style(
      HTML(".shiny-notification{position:fixed;top:calc(2.5%);left:calc(60%);
           width:280px;height:50px;font-size:17px;opacity:1;font-family:'Verdana',sans-serif}"))),
  tags$head(tags$style(".modal-dialog.xl1{width:95%}")),
  tags$head(tags$meta(charset="UTF-8")),

  # titlePanel("SCIP: Suite for CNV Interpretation and Prioritization"),
  
  fluidRow(
    column(6,selectInput(inputId="cnv_name",label="CNV名称 (开始键入以搜索)", 
              choices=NULL,width="100%")),
    column(2,actionButton(inputId="save",label="保存",width="100%"),style="margin-top:25px;"),
    column(2,actionButton(inputId="prev_var",label="上一个",width="100%"),style="margin-top:25px;"),
    column(2,actionButton(inputId="next_var",label="下一个",width="100%"),style="margin-top:25px"),
  ),
  
  fluidRow(
    column(2,selectInput(inputId="qual_override",label="变体质量",
           choices=c("无覆盖"="No Override","通过"="Passed","失败"="Failed"),selected="No Override",width="100%")),
    
    column(4,selectInput(inputId="decision",label="测定",
                         choices=c("未评估"="Not Evaluated","排除-质量不足/难以评估"="Ruled Out - Quality Inadequate / Difficult to Assess","排除-人口变异"="Ruled Out - Population Variation",
                                   "排除-没有发现感兴趣的基因"= "Ruled Out - No Gene of Interest Identified","排除-不正确的边界，完全内含子"="Ruled Out - Incorrect Boundary, Fully Intronic",
                                   "排除-非基因内DUP"="Ruled Out - Non-intragenic DUP","排除-其他原因"="Ruled Out - Other Reasons",
                                   "延迟-隐性基因-寻找复合Het SNV"= "Deferred - Recessive Gene - Look for Compound Het SNV", "已解释-由另一个调用方标识的相同变体"="Already Interpreted - Same Variant Identified by Another Caller",
                                   "进一步审查-可能需要报告"="Further Review - Potentially Reportable","进一步审查-不太可能报告"="Further Review - Not Likely Reportable"),selected="Not Evaluated",width="100%")),
    
    column(6,textAreaInput(inputId="comment",label="注释 (可选)",placeholder="关于变体的可选说明",
                       width="100%")),
  ),
  textOutput("last_edit_info"),
  hr(),
  
  h4("变体摘要"), 
  DT::dataTableOutput("basic_info",width="1250px"),
  DT::dataTableOutput("basic_links",width="250px"),
  br(),
  DT::dataTableOutput("highlights1",width="800px"),
  DT::dataTableOutput("highlights2",width="800px"),
  
  hr(),
  h4("读取深度和映射质量"),
  fluidRow(
    column(2,selectInput(inputId="binsize",label="bin大小",choices=c("1 bp"=1,"10 bp"=10,"25 bp"=25,"50 bp"=50,"100 bp"=100,"250 bp"=250,"500 bp"=500,
                                                                      "1 kb"=1000,"5 kb"=5000,"10 kb"=10000,"50 kb"=50000,"100 kb"=99999.999),selected="1000",width="100%")),
    column(2,actionButton(inputId="render_depth",label="更新读取深度",width="100%"),style="margin-top:25px;"),
    column(3,selectInput(inputId="use_12878",label="NA12878",choices=c("仅查看示例"="sample_only","查看覆盖了NA12878的示例"="both","仅查看NA12878"="12878_only"),
                         selected="both",width="100%")),
    column(2,actionButton(inputId="render_mq",label="更新映射质量",width="100%"),style="margin-top:25px;"),
  ),
  
  plotOutput('depth',dblclick='click_depth',brush=brushOpts(id="brush_depth",resetOnNew=T),height="250px"),
  plotOutput('mq',dblclick='click_mq',brush=brushOpts(id="brush_mq",resetOnNew=T),height="250px"),
  
  hr(),
  h4("异常读取"),
  fluidRow(
    column(2,actionButton(inputId="render_reads",label="加载异常读取",width="100%"),style="margin-top:25px;"),
    column(3,selectInput(inputId="read_type",label="读取类型",choices=c("全部 (双端和拆分读取)"="all","双端只读"="pe",
                                                                    "正常方向，仅小插入尺寸"="si","正常方向，仅大插入尺寸"="pe_abnormal_orientation",
                                                                    "拆分-仅读取"="sr"),selected="all",width="100%")),
    column(3,sliderInput(inputId="pe_percentile",label="双端异常值百分位数",value=99.5,min=90,max=99.9,step=0.1)),
  ),
  p("插入尺寸估计值"),
  tableOutput("insert_size_estimates"),
  plotOutput("reads",dblclick='click_reads',brush=brushOpts(id="brush_reads",resetOnNew=T),height="400px"),
  DT::dataTableOutput("reads_table",width="100%"),
  
  hr(),
  h4("外部和内部变体数据库"),
  fluidRow(
    column(3,selectInput(inputId="gnomad_af",label="gnomAD SV等位基因频率",choices=c("No Filter"=0,"Above 0.01%"=1e-4,"Above 0.05%"=5e-4,"Above 0.1%"=1e-3,
                                                                                                                "Above 0.5%"=5e-3,"Above 1%"=0.01,"Above 5%"=0.05),selected=0,width="100%")),
    column(2,selectInput(inputId="clinvar_subset",label="Clinvar-后果",choices=c("无筛选器"="all","仅P/LP"="p","任何非B/LB"="non-b",
                                                                                          "仅B/LB"="b"),selected="all",width="100%")),
    column(2,selectInput(inputId="clinvar_size",label="Clinvar-变体尺寸",choices=c("No Filter"=1e20,"10 Mb以下"=10e6,"5 Mb以下"=5e6,"1 Mb以下"=1e6,
                                                                                         "500 Kb以下"=5e5,"250 Kb以下"=2.5e5),selected=5e6,width="100%")),
  ),
  plotOutput("gnomadsv_clinvar",dblclick='click_gnomadsv_clinvar',brush=brushOpts(id="brush_gnomadsv_clinvar",resetOnNew=T),height="600px"),
  br(),
  tabsetPanel(type="tabs",
              tabPanel("gnomAD SV",DT::dataTableOutput("gnomadsv_table",width="800px")),
              tabPanel("ClinVar",DT::dataTableOutput("clinvar_table",width="100%")),
              tabPanel("内部队列",DT::dataTableOutput("cgc_table",width="800px")),
              selected="gnomAD SV"
  ),
  
  hr(),
  h4("基因组邻域"),
  fluidRow(
    column(3,actionButton(inputId="render_genes",label="绘制基因组邻域",width="100%"))),
  plotOutput("genes",dblclick='click_genes',brush=brushOpts(id="brush_genes",resetOnNew=T),height="400px"),
  tabsetPanel(type="tabs",
              tabPanel("基因",DT::dataTableOutput("genes_table",width="100%")),
              tabPanel("Clingen剂量图",DT::dataTableOutput("dosage_table",width="1000px")),
              selected="基因"
  ),

  hr(),
  h4("重要说明"),
  p("对于Manta INS和BND变体: 读取深度和映射质量部分不适用。异常读取部分是实验性的。鼓励用户在IGV中仔细检查这些变体。"),
  p("在外部和内部变体数据库部分中，未显示在所分析的变体的相反方向上的CNV (例如，当所分析的变体是DUP时的DEL)。"),
  p("在基因组邻域部分中，未显示与所分析的变体不重叠的Clingen HI/TS区域/基因。对于DELS，不显示TS信息。"),
    
  hr(),
  h4("版本控制"),
  # tableOutput("version_control1"),
  tableOutput("version_control2"),
  
  # hr(),
  # p("Interface Version 0.1.7 Public (20220202)"),
  # p("Cardiac Genome Clinic, Ted Rogers Centre for Heart Research"),
  # p("Division of Clinical and Metabolic Genetics & The Centre for Applied Genomics, The Hospital for Sick Children. © 2022"),
)

get_scip_cnv_name <- function(sample_id, cnv_param) {
  # Split the input by "-"
  parts <- unlist(strsplit(cnv_param, "-"))
  
  # Extract components
  chrom <- parts[1]     # 13
  start <- parts[2]     # 35698501
  end <- parts[3]       # 35731500
  cnv_type <- parts[4]  # DEL
  
  # Construct the new format
  cnv_name <- sprintf("%s-001-001.%s.%s.%s.%s", sample_id, chrom, start, end, cnv_type)
  
  return(cnv_name)
}

server=function(input,output,session){
  # Reactive values to store dynamic data
  rv <- reactiveValues(
    config_file = NULL,
    temp_file_dir = NULL,
    username = NULL,
    regions_list = NULL,
    regions_name = NULL,
    config_loaded = FALSE,
    cnv_name = NULL,
    sample_id = NULL,
    assembly = NULL
  )

  # Read URL parameters and update config
  observe({
    query <- parseQueryString(session$clientData$url_search)
    
    # Use config file from URL, or fallback to default
    # Default config path
    config_path <- "./interface_config.txt"

    # If sample_id exists in the URL, update config_path and store sample_id
    if (!is.null(query$sample_id)) {
      rv$sample_id <- query$sample_id
      rv$results_out <- paste0("/glsofort/samples/", query$sample_id, "/SCIP_results.txt")
      config_path <- paste0("/highspeed-data/samples/", query$sample_id, "/SCIP/interface_config.txt")
    }
    
    rv$assembly <- "hg19"
    if (!is.null(query$sample_id)) {
      rv$assembly <- query$assembly
    }

    # Read config file dynamically
    rv$config_file <- read.table(config_path, comment.char = "", sep = "\t", quote = "")
    
    # Extract file paths from config
    rv$temp_file_dir <- rv$config_file[2, 2]
    rv$username <- rv$config_file[4, 2]
    
    # Load CNV regions list dynamically
    rv$regions_list <- read.table(
      paste0(rv$temp_file_dir, "/../user_data/", rv$config_file[1, 2], ".pipeline_summary.txt"),
      comment.char = "", sep = "\t", quote = ""
    )

    # Sort and update regions_name dynamically
    regions_sorted <- rv$regions_list[sort(rv$regions_list$V5 - rv$regions_list$V4 + 1, 
                                           index.return = TRUE, decreasing = TRUE)$ix,]
    rv$regions_name <- regions_sorted$V17[sort(regions_sorted$V18, index.return = TRUE)$ix]
    if (file.exists(rv$results_out)==F){
      file.create(rv$results_out)
    }

    # Get CNV
    if (is.null(query$cnv)) {
      rv$cnv_name <- rv$regions_name[1]
    } else {
      print(rv$sample_id)
      rv$cnv_name <- get_scip_cnv_name(rv$sample_id, query$cnv)
      print(rv$cnv_name)
    }

    updateSelectInput(session, "cnv_name", choices = rv$regions_name)
    updateSelectInput(session, "cnv_name", choices = rv$regions_name, selected = rv$cnv_name) 
    print(rv$cnv_name)

    # Mark config as ready
    rv$config_loaded <- TRUE

  })

  # reactive values
  ranges=reactiveValues(x=NULL,y1=NULL,name=NULL)
  ranges_default=reactiveValues(x=NULL,y1=NULL)
  load=reactiveValues(depth=0,mq=0,reads=0,gene=0)
  master_info=reactiveValues(chr=NULL,start=NULL,end=NULL,proband=NULL,type=NULL,loaded=FALSE)
  
  latest_interpretation=reactiveValues(qual_override=NULL,decision=NULL,comment=NULL,time="从不",user="无")
  insert_size=reactiveValues(lower=NULL,upper=NULL,name=NULL,percentile=NULL)
  abnormal_reads=reactiveValues(pe_orientation=NULL,pe_small=NULL,pe_large=NULL,se=NULL)
  out_reads=reactiveValues(table=NULL)
  sv1=reactiveValues(gnomad=NULL,clinvar=NULL,cgc=NULL,dosage_hi=NULL,dosage_ts=NULL)
  # end reactive values
  
  observeEvent(input$prev_var,{
    req(rv$config_loaded)
    current_var=which(rv$regions_name==input$cnv_name)
    st2=0
    if (current_var>1){
      updateSelectInput(inputId="cnv_name",selected=rv$regions_name[current_var-1])
    } else{
      showNotification("已经在第一个变体",type="error",duration=15)
      st2=1
    }
    
    if ((input$qual_override!="No Override" && input$qual_override!=latest_interpretation$qual_override) || 
        (input$decision!="Not Evaluated" && input$decision!=latest_interpretation$decision) || (input$comment!="" && input$comment!=latest_interpretation$comment)){
      write(paste(as.numeric(Sys.time()),input$cnv_name,input$qual_override,input$decision,input$comment,Sys.time(),rv$username,sep="\t"),
            file=rv$results_out,append=T)
      if (st2==0){
        showNotification("已保存解释",type="message",duration=2)
      } else{
        latest_interpretation$qual_override=input$qual_override
        latest_interpretation$decision=input$decision
        latest_interpretation$comment=input$comment
        latest_interpretation$time=Sys.time()
        latest_interpretation$user=rv$username
        showNotification("已保存解释",type="message",duration=1)
      }
    }
  })
  
  observeEvent(input$next_var,{
    req(rv$config_loaded)
    current_var=which(rv$regions_name==input$cnv_name)
    st2=0
    if (current_var<length(rv$regions_name)){
      updateSelectInput(inputId="cnv_name",selected=rv$regions_name[current_var+1])
    } else{
      showNotification("已经在最后一个变体",type="error",duration=15)
      st2=1
    }
    
    if ((input$qual_override!="No Override" && input$qual_override!=latest_interpretation$qual_override) || 
        (input$decision!="Not Evaluated" && input$decision!=latest_interpretation$decision) || (input$comment!="" && input$comment!=latest_interpretation$comment)){
      write(paste(as.numeric(Sys.time()),input$cnv_name,input$qual_override,input$decision,input$comment,Sys.time(),rv$username,sep="\t"),
            file=rv$results_out,append=T)
      if (st2==0){
        showNotification("已保存解释",type="message",duration=2)
      } else{
        latest_interpretation$qual_override=input$qual_override
        latest_interpretation$decision=input$decision
        latest_interpretation$comment=input$comment
        latest_interpretation$time=Sys.time()
        latest_interpretation$user=rv$username
        showNotification("已保存解释",type="message",duration=1)
      }
    }
  })
  
  observeEvent(input$save,{
    req(rv$config_loaded)
    if ((input$qual_override!="No Override" && input$qual_override!=latest_interpretation$qual_override) || 
        (input$decision!="Not Evaluated" && input$decision!=latest_interpretation$decision) || (input$comment!="" && input$comment!=latest_interpretation$comment)){
      write(paste(as.numeric(Sys.time()),input$cnv_name,input$qual_override,input$decision,input$comment,Sys.time(),rv$username,sep="\t"),
            file=rv$results_out,append=T)
      latest_interpretation$qual_override=input$qual_override
      latest_interpretation$decision=input$decision
      latest_interpretation$comment=input$comment
      latest_interpretation$time=Sys.time()
      latest_interpretation$user=rv$username
      showNotification("已保存解释",type="message",duration=2)
    } else{
      showNotification("与上次保存时无更改",type="message",duration=2)
    }
  })

  observeEvent(input$cnv_name,{
    req(rv$config_loaded, input$cnv_name)
    info=as.data.frame(strsplit(input$cnv_name,"\\."))[,1]
    master_info$proband=info[1]
    master_info$chr=info[2]
    master_info$start=as.numeric(info[3])
    master_info$end=as.numeric(info[4])
    master_info$type=info[5]
    master_info$loaded <- TRUE

    updateSelectInput(inputId="use_12878",selected="both")
    updateSelectInput(inputId="read_type",selected="all")
    updateSliderInput(inputId="pe_percentile",value=99.5)
    updateSelectInput(inputId="gnomad_af",selected=0)
    updateSelectInput(inputId="clinvar_subset",selected="all")
    updateSelectInput(inputId="clinvar_size",selected=5e6)
    
    load$depth=0
    load$mq=0
    load$reads=0
    load$genes=0
    
    length=as.numeric(master_info$end)-as.numeric(master_info$start)+1
    if (length<=10000){
      updateSelectInput(inputId="binsize",selected="250")
    } else if (length<=50000){
      updateSelectInput(inputId="binsize",selected="500")
    } else if (length<=500000){
      updateSelectInput(inputId="binsize",selected="1000")
    } else if (length<=1000000){
      updateSelectInput(inputId="binsize",selected="5000")
    } else if (length<=2000000){
      updateSelectInput(inputId="binsize",selected="10000")
    } else if (length<=5000000){
      updateSelectInput(inputId="binsize",selected="50000")
    } else{
      updateSelectInput(inputId="binsize",selected="99999.999")
    }
    
    current_var=which(rv$regions_name==input$cnv_name)
    x1=try(read.table(rv$results_out,comment.char="",sep="\t",quote=""),silent=T)
    st1=0
    if(class(x1)=="data.frame"){
      x1_sub=x1[which(x1$V2==rv$regions_name[current_var]),]
      if(length(x1_sub[,1])>0){
        x1_sub_sorted=x1_sub[sort(x1_sub$V1,index.return=T,decreasing=T)$ix,]
        updateSelectInput(inputId="qual_override",selected=x1_sub_sorted[1,3])
        updateSelectInput(inputId="decision",selected=x1_sub_sorted[1,4])
        updateTextAreaInput(inputId="comment",value=x1_sub_sorted[1,5])
        latest_interpretation$qual_override=x1_sub_sorted[1,3]
        latest_interpretation$decision=x1_sub_sorted[1,4]
        latest_interpretation$comment=x1_sub_sorted[1,5]
        latest_interpretation$time=x1_sub_sorted[1,6]
        latest_interpretation$user=x1_sub_sorted[1,7]
        st1=1
      }
    }
    if(st1==0){
      updateSelectInput(inputId="qual_override",selected="No Override")
      updateSelectInput(inputId="decision",selected="Not Evaluated")
      updateTextAreaInput(inputId="comment",value="")
      latest_interpretation$qual_override="No Override"
      latest_interpretation$decision="Not Evaluated"
      latest_interpretation$comment=""
      latest_interpretation$time="从不"
      latest_interpretation$user="无"
    }
    
    if(is.null(ranges$x)==T || ranges$name!=input$cnv_name){
      x2=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script05_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
      ranges$name=input$cnv_name
      ranges$x=c(min(x2$V1,na.rm=T),max(x2$V1,na.rm=T))
      ranges$y1=c(0,quantile(x2$V2,.95,na.rm=T)*1.5)
      ranges_default$x=c(min(x2$V1,na.rm=T),max(x2$V1,na.rm=T))
      ranges_default$y1=c(0,quantile(x2$V2,.95,na.rm=T)*1.5)
    }
  })
  
  # action buttons for loading specific section of the webpage
  observeEvent(input$render_depth,{
    req(rv$config_loaded)
    load$depth=load$depth+1
  })
  
  observeEvent(input$render_mq,{
    req(rv$config_loaded)
    load$mq=load$mq+1
  })
  
  observeEvent(input$render_reads,{
    req(rv$config_loaded)
    load$reads=load$reads+1
  })
  
  observeEvent(input$render_genes,{
    req(rv$config_loaded)
    load$genes=load$genes+1
  })
  ### end action buttons
  
  output$last_edit_info=renderText(paste("最后解释: ",latest_interpretation$time," (",latest_interpretation$user,")",sep=""))
  
  output$basic_info=DT::renderDataTable({
    req(rv$config_loaded, input$cnv_name)
    x3=rv$regions_list[which(rv$regions_list$V17==input$cnv_name),]
    auto_result=as.data.frame(strsplit(as.character(x3[1]),"\\s"))[1,1]
    length=x3[5]-x3[4]+1

    if (x3[18]<99){
      priority_print=paste(x3[18],"  ","<span class='label label-danger'>高</span>",sep="")
    }
    else if (x3[18]==99){
      priority_print=paste(x3[18],"  ","<span class='label label-warning'>中</span>",sep="")
    }
    else{
      priority_print=paste(x3[18],"  ","<span class='label label-info'>低</span>",sep="")
    }
    
    as.data.frame(rbind(c("质量","优先级","ID","Chr","开始","结束","类型","大小 (kb)","深度比","映射质量",
                    "支撑对","对立对","拆分-读取"),
                  c(auto_result,priority_print,x3[2],x3[3],x3[4],x3[5],x3[6],round(length*1e3/1e3)/1e3,round(x3[9]*1e3)/1e3,round(x3[10]*1e3)/1e3,x3[13],x3[14],x3[16])))
  },options=list(dom="t",language=list(info="显示 _START_ 个条目中的 _END_ 到 _TOTAL_ 个",paginate=list(previous="上一个",`next`="下一个"),search = "搜索:",lengthMenu ="显示 _MENU_ 个条目"),ordering=F,headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F,selection='none')
  
  output$basic_links=DT::renderDataTable({
    req(rv$config_loaded, master_info$loaded)
    if (rv$assembly == "hg19") {
      as.data.frame(cbind(paste("<a href='https://www.deciphergenomics.org/browser#q/grch37:",master_info$chr,":",master_info$start,"-",master_info$end,"/location/grch37' target='_blank'>DECIPHER</a>",sep=""),
                    paste("<a href='https://gnomad.broadinstitute.org/region/",master_info$chr,"-",master_info$start,"-",master_info$end,"?dataset=gnomad_sv_r2_1' target='_blank'>gnomAD SV</a>",sep=""),
                    paste("<a href='http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19/?name=chr",master_info$chr,"%3A",master_info$start,"-",master_info$end,";search=Search' target='_blank'>DGV</a>",sep="")))
    } else {
      as.data.frame(cbind(paste("<a href='https://www.deciphergenomics.org/browser#q/",master_info$chr,":",master_info$start,"-",master_info$end,"/location/grch38' target='_blank'>DECIPHER</a>",sep=""),
                    paste("<a href='http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg38/?name=chr",master_info$chr,"%3A",master_info$start,"-",master_info$end,";search=Search' target='_blank'>DGV</a>",sep="")))
    }
  },options=list(dom="t",language=list(info="显示 _START_ 个条目中的 _END_ 到 _TOTAL_ 个",paginate=list(previous="上一个",`next`="下一个"),search = "搜索:",lengthMenu ="显示 _MENU_ 个条目"),ordering=F,headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F,escape=F,selection='none')
  
  output$highlights1=DT::renderDataTable({
    req(rv$config_loaded, input$cnv_name, master_info$loaded)
    x18=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script08_file2.txt.gz",sep=""),comment.char="",sep="\t",quote="")
    if (length(x18$V1[which(x18$V2==1)])>=1){
      datatable(rbind("正面信息:",as.data.frame(x18$V1[which(x18$V2==1)])),class="compact",selection='none',
                options=list(dom="t",language=list(info="显示 _START_ 个条目中的 _END_ 到 _TOTAL_ 个",paginate=list(previous="上一个",`next`="下一个"),search = "搜索:",lengthMenu ="显示 _MENU_ 个条目"),ordering=F,headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F,escape=F) %>% DT::formatStyle(1,color="#ef6548",fontSize="110%")
    } else{
      datatable(as.data.frame("无正面信息"),class="compact",selection='none',
                options=list(dom="t",ordering=F,headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F) %>% DT::formatStyle(1,color="#ef6548",fontSize="110%")
    }
  }) 
  
  output$highlights2=DT::renderDataTable({
    req(rv$config_loaded, input$cnv_name, master_info$loaded)
    x18=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script08_file2.txt.gz",sep=""),comment.char="",sep="\t",quote="")
    if (length(x18$V1[which(x18$V2==2)])>=1){
      datatable(rbind("负面信息:",as.data.frame(x18$V1[which(x18$V2==2)])),class="compact",selection='none',
                options=list(dom="t",language=list(info="显示 _START_ 个条目中的 _END_ 到 _TOTAL_ 个",paginate=list(previous="上一个",`next`="下一个"),search = "搜索:",lengthMenu ="显示 _MENU_ 个条目"),ordering=F,headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F) %>% DT::formatStyle(1,color="#02818a",fontSize="110%")
    } else{
      datatable(as.data.frame("无负面信息"),class="compact",selection='none',
                options=list(dom="t",ordering=F,headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F) %>% DT::formatStyle(1,color="#02818a",fontSize="110%")  
    }
  })
  
  output$depth=renderPlot({
    req(rv$config_loaded, input$cnv_name, master_info$loaded)
    par(mgp=c(2.1,1,0),mai=c(0.4,0.62,0.5,0.1))
    length=as.numeric(master_info$end)-as.numeric(master_info$start)+1
    x6=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script05_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
    
    plot(0,0,col=rgb(0,0,0,0),
         xlab="",xlim=ranges$x,ylab="读取深度",ylim=ranges$y1,xaxs="i",yaxs="i",main="读取深度")
    if (load$depth==0){
      par(new=T)
      plot(0,0,col=rgb(0,0,0,0),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
      text(0.5,0.5,"单击上面的 \"更新读取深度\" 按钮以显示读取深度。\n
           对于较大的变体，建议将bin大小设置为5 kb或更大 (在单击按钮之前)，以提高性能。",cex=1.5)
    } else if (master_info$type=="BND" || master_info$type=="INS"){
      plot(0,0,col=rgb(0,0,0,0),
           xlab="",xlim=ranges$x,ylab="读取深度",ylim=ranges$y1,xaxs="i",yaxs="i",main="读取深度")
      par(new=T)
      plot(0,0,col=rgb(0,0,0,0),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
      text(0.5,0.5,"本节不适用于INS或BND类型变体。",cex=1.5)
      box()
    } else if (load$depth>0){
      # speed up using BinMean and plotH
      # using isolate so the plot will need to be manually updated
      x6_sub=x6[intersect(which(x6$V1>=ranges$x[1]-as.numeric(isolate(input$binsize))*2),which(x6$V1<=ranges$x[2]+as.numeric(isolate(input$binsize))*2)),]
      
      if (input$use_12878=="both" || input$use_12878=="12878_only"){
        x7=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script04_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
        norm_rd_12878=median(x7$V2[c(which(x7$V1<=master_info$start),which(x7$V1>=master_info$end))],na.rm=T)/median(x6$V2[c(which(x6$V1<=master_info$start),which(x6$V1>=master_info$end))],na.rm=T)
        x7_sub=x7[intersect(which(x7$V1>=ranges$x[1]-as.numeric(isolate(input$binsize))*2),which(x7$V1<=ranges$x[2]+as.numeric(isolate(input$binsize))*2)),]
        par(new=T)
        plot(0,0,col=rgb(0,0,0,0),
             xlab="",ylab="",xlim=ranges$x,ylim=ranges$y1,xaxs="i",yaxs="i",axes=F)
        
        # speed up using BinMean and plotH
        x7_pos=BinMean(x7_sub$V1,every=as.numeric(isolate(input$binsize)))
        x7_val=BinMean(x7_sub$V2,every=as.numeric(isolate(input$binsize)))/norm_rd_12878
        median_out_x7=median(x7$V2[c(which(x7$V1<=master_info$start),which(x7$V1>=master_info$end))],na.rm=T)/norm_rd_12878
        median_in_x7=median(x7$V2[intersect(which(x7$V1>master_info$start),which(x7$V1<master_info$end))],na.rm=T)/norm_rd_12878
        par(new=T)
        plotH(x7_pos,x7_val,width=as.numeric(isolate(input$binsize)),col="#ff7f0080",xlab="",ylab="",xlim=ranges$x,ylim=ranges$y1,xaxs="i",yaxs="i",axes=F,border=NA)
        lines(x=c(ranges_default$x[1],master_info$start),y=c(median_out_x7,median_out_x7),col="#d95f0e",lty=2)
        lines(x=c(master_info$end,ranges_default$x[2]),y=c(median_out_x7,median_out_x7),col="#d95f0e",lty=2)
        lines(x=c(master_info$start,master_info$end),y=c(median_in_x7,median_in_x7),col="#d95f0e",lty=2)        
        # for(i in seq(1,length(x7$V1),as.numeric(input$binsize))){
        #   end=i+as.numeric(input$binsize)-1
        #   rect(x7$V1[i]-0.5,0,x7$V1[end]+0.5,mean(x7$V2[i:end],na.rm=T)/norm_rd_12878,border=NA,col="#ff7f0080")
        # }
        if (input$use_12878=="both"){
          legend("topleft",legend=c(master_info$proband,"NA12878"),fill=c("#96969680","#ff7f0080"),cex=0.8,bg="white",border=NA)
        }
      }
      
      if(input$use_12878=="both" || input$use_12878=="sample_only"){
        x6_pos=BinMean(x6_sub$V1,every=as.numeric(isolate(input$binsize)))
        x6_val=BinMean(x6_sub$V2,every=as.numeric(isolate(input$binsize)))
        median_out_x6=median(x6$V2[c(which(x6$V1<=master_info$start),which(x6$V1>=master_info$end))],na.rm=T)
        median_in_x6=median(x6$V2[intersect(which(x6$V1>master_info$start),which(x6$V1<master_info$end))],na.rm=T)
        par(new=T)
        plotH(x6_pos,x6_val,width=as.numeric(isolate(input$binsize)),col="#96969680",xlab="",ylab="",xlim=ranges$x,ylim=ranges$y1,xaxs="i",yaxs="i",axes=F,border=NA)
        lines(x=c(ranges_default$x[1],master_info$start),y=c(median_out_x6,median_out_x6),col="#636363",lty=2)
        lines(x=c(master_info$end,ranges_default$x[2]),y=c(median_out_x6,median_out_x6),col="#636363",lty=2)
        lines(x=c(master_info$start,master_info$end),y=c(median_in_x6,median_in_x6),col="#636363",lty=2)
      }
      
      # for(i in seq(1,length(x6$V1),as.numeric(input$binsize))){
      #   end=i+as.numeric(input$binsize)-1
      #   rect(x6$V1[i]-0.5,0,x6$V1[end]+0.5,mean(x6$V2[i:end],na.rm=T),border=NA,col="#969696")
      # }
      
      rect(master_info$start,-1000,master_info$end,1000,col="#8dd3c733",border="#8dd3c74D")
    }
  })
  
  output$mq=renderPlot({
    req(rv$config_loaded, input$cnv_name, master_info$loaded)
    par(mgp=c(2.1,1,0),mai=c(0.4,0.62,0.5,0.1))

    plot(0,0,col=rgb(0,0,0,0),
         xlab="",xlim=ranges$x,ylab="映射质量",ylim=c(0,65),xaxs="i",yaxs="i",main="映射质量")
    if (load$mq==0){
      par(new=T)
      plot(0,0,col=rgb(0,0,0,0),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
      text(0.5,0.5,"单击上面的 \"更新映射质量\" 按钮以显示映射质量。\n
           对于较大的变体，建议将bin大小设置为5 kb或更大 (在单击按钮之前)，以提高性能。",cex=1.5)
    } else if (master_info$type=="BND" || master_info$type=="INS"){
      plot(0,0,col=rgb(0,0,0,0),
           xlab="",xlim=ranges$x,ylab="映射质量",ylim=c(0,65),xaxs="i",yaxs="i",main="映射质量")
      par(new=T)
      plot(0,0,col=rgb(0,0,0,0),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
      text(0.5,0.5,"本节不适用于INS或BND类型变体。",cex=1.5)
      box()
    } else if (load$mq>0){
      if (input$use_12878=="both" || input$use_12878=="12878_only"){
        x7=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script04_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
        x7_sub=x7[intersect(which(x7$V1>=ranges$x[1]-as.numeric(isolate(input$binsize))*2),which(x7$V1<=ranges$x[2]+as.numeric(isolate(input$binsize))*2)),]
        
        # speed up using BinMean and plotH
        x7_pos=BinMean(x7_sub$V1,every=as.numeric(isolate(input$binsize)))
        x7_val=BinMean(x7_sub$V3,every=as.numeric(isolate(input$binsize)))
        par(new=T)
        plotH(x7_pos,x7_val,width=as.numeric(isolate(input$binsize)),col="#ff7f0080",xlab="",ylab="",xlim=ranges$x,ylim=c(0,65),xaxs="i",yaxs="i",axes=F,border=NA)
        if (input$use_12878=="both"){
          legend("topleft",legend=c(master_info$proband,"NA12878"),fill=c("#96969680","#ff7f0080"),cex=0.8,bg="white",border=NA)
        }
      }
      
      if(input$use_12878=="both" || input$use_12878=="sample_only"){
        x6=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script05_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
        x6_sub=x6[intersect(which(x6$V1>=ranges$x[1]-as.numeric(isolate(input$binsize))*2),which(x6$V1<=ranges$x[2]+as.numeric(isolate(input$binsize))*2)),]
        
        # speed up using BinMean and plotH
        x6_pos=BinMean(x6_sub$V1,every=as.numeric(isolate(input$binsize)))
        x6_val=BinMean(x6_sub$V3,every=as.numeric(isolate(input$binsize)))
        par(new=T)
        plotH(x6_pos,x6_val,width=as.numeric(isolate(input$binsize)),col="#96969680",xlab="",ylab="",xlim=ranges$x,ylim=c(0,65),xaxs="i",yaxs="i",axes=F,border=NA)
      }
      rect(master_info$start,-1000,master_info$end,1000,col="#8dd3c733",border="#8dd3c74D")
    }
  })
  
  observeEvent(input$click_depth,{
    req(rv$config_loaded)
    brush=input$brush_depth
    if (!is.null(brush)){
      ranges$x=c(brush$xmin,brush$xmax)
      ranges$y1=c(brush$ymin,brush$ymax)
    }else{
      ranges$x=ranges_default$x
      ranges$y1=ranges_default$y1
    }
  })
  
  observeEvent(input$click_mq,{
    req(rv$config_loaded)
    brush=input$brush_mq
    if (!is.null(brush)){
      ranges$x=c(brush$xmin,brush$xmax)
    }else{
      ranges$x=ranges_default$x
    }
  })
  
  output$insert_size_estimates=renderTable({
    req(rv$config_loaded, input$cnv_name)
    if (load$reads>0){
      x8=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script05_file2.txt.gz",sep=""),comment.char="",sep="\t",quote="")
      # added to remove read pairs/split-reads that have one part completely outside the extended CNV region
      # x8=x8[intersect(which(x8$V5>=ranges_default$x[1]),which(x8$V6<=ranges_default$x[2])),]
      
      all_size=x8$V8[which(x8$V11=="FR_inward")]
      if (is.null(insert_size$lower)==T || insert_size$name!=input$cnv_name || insert_size$percentile!=input$pe_percentile){
        insert_size$lower=quantile(all_size,1-input$pe_percentile/100)
        insert_size$higher=quantile(all_size,input$pe_percentile/100)
        insert_size$name=input$cnv_name
        insert_size$percentile=input$pe_percentile
      }
      
      abnormal_reads$pe_orientation=x8[intersect(which(x8$V11!="FR_inward"),which(x8$V11!="SR")),]
      abnormal_reads$pe_large=x8[intersect(which(x8$V11=="FR_inward"),which(x8$V8>=insert_size$higher)),]
      abnormal_reads$pe_small=x8[intersect(which(x8$V11=="FR_inward"),which(x8$V8<=insert_size$lower)),]
      abnormal_reads$sr=x8[which(x8$V11=="SR"),]
      as.data.frame(rbind(c("下界",round(insert_size$lower*1e2)/1e2),c("上界",round(insert_size$higher*1e2)/1e2)))
    } else{
      as.data.frame(c("单击上面的 \"加载异常读数\" 按钮以加载此部分。"))
    }
  },colnames=F,spacing="s",digits=2,striped=T,bordered=T,align="c")
  
  output$reads=renderPlot({
    req(rv$config_loaded, master_info$loaded, input$read_type)
    par(mgp=c(2.1,1,0),mai=c(0.4,0.62,0.5,0.1))
    if (load$reads==0){
      plot(0,0,col=rgb(0,0,0,0),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="",main="异常读取")
      box()
      text(0.5,0.5,"单击上面的 \"加载异常读数\" 按钮以加载此部分。\n\nINS和BND型变体的实验。",cex=2)
    }
    
    if (load$reads>0){
      table_out=NULL
      if (input$read_type=="all"){
        table_out=rbind(abnormal_reads$pe_orientation,abnormal_reads$pe_small,abnormal_reads$pe_large,abnormal_reads$sr)
      } else if (input$read_type=="pe"){
        table_out=rbind(abnormal_reads$pe_orientation,abnormal_reads$pe_small,abnormal_reads$pe_large)
      } else if (input$read_type=="pe_abnormal_orientation"){
        table_out=abnormal_reads$pe_orientation
      } else if (input$read_type=="si"){
        table_out=abnormal_reads$pe_small
      } else if (input$read_type=="li"){
        table_out=abnormal_reads$pe_large
      } else if (input$read_type=="sr"){
        table_out=abnormal_reads$sr
      }
      colnames(table_out)=c("Read_Name","Flag","Chr","Left_Start","Left_End","Right_Start","Right_End","Size","MQ","CIGAR","Read_Type")

      # table_out_interval=table_out[union(intersect(which(table_out$Left_Start>=ranges$x[1]),which(table_out$Right_End<=ranges$x[2])),
      #                                    intersect(which(table_out$Left_Start>=ranges$x[1]),which(table_out$Right_End<=ranges$x[2]))),]
      # table_out_interval=table_out[intersect(which(table_out$Left_Start<=ranges$x[2]),which(table_out$Right_End>=ranges$x[1])),]
      
      table_out_interval=table_out[union(intersect(which(table_out$Left_End>=ranges$x[1]),which(table_out$Left_Start<=ranges$x[2])),
                                         intersect(which(table_out$Right_End>=ranges$x[1]),which(table_out$Right_Start<=ranges$x[2]))),]
      table_out_interval=table_out_interval[sort(table_out_interval$Left_Start,index.return=T,decreasing=T)$ix,]
      table_out_interval$Read_Type[which(table_out_interval$Read_Type=="SR")]="Split-read"
      out_reads$table=table_out_interval
      
      if (length(table_out_interval$Chr)>=1){
        plot(0,0,col=rgb(0,0,0,0),xlim=ranges$x,ylim=c(0,(length(table_out_interval$Chr)+1)),axes=F,xlab="",ylab="",xaxs="i",yaxs="i",main="异常读取")
        box()
        axis(1)
  
        for (i in 1:length(table_out_interval$Chr)){
          if (length(table_out_interval$Chr)<=50){
            read_name=unlist(strsplit(table_out_interval$Read_Name[i],"\\:"))
            text(table_out_interval$Right_End[i],i+0.1,read_name[length(read_name)],pos=4,cex=0.8)
          }
          
          if(table_out_interval$Read_Type[i]=="FR_inward"){
            if (table_out_interval$Left_End[i]-25>=table_out_interval$Left_Start[i]){
              rect(table_out_interval$Left_Start[i],i-0.25,table_out_interval$Left_End[i]-25,i+0.25,border=NA,col="#e41a1cBF")
            }
            if (table_out_interval$Right_Start[i]+25<=table_out_interval$Right_End[i]){
              rect(table_out_interval$Right_Start[i]+25,i-0.25,table_out_interval$Right_End[i],i+0.25,border=NA,col="#e41a1cBF")
            }
            polygon(x=c(table_out_interval$Left_End[i]-25,table_out_interval$Left_End[i],table_out_interval$Left_End[i]-25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#e41a1cBF")
            polygon(x=c(table_out_interval$Right_Start[i]+25,table_out_interval$Right_Start[i],table_out_interval$Right_Start[i]+25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#e41a1cBF")
            lines(x=c(table_out_interval$Left_End[i],table_out_interval$Right_Start[i]),y=c(i,i),col="#e41a1cBF")
  
          } else if (table_out_interval$Read_Type[i]=="FR_outward"){
            if (table_out_interval$Left_Start[i]+25<=table_out_interval$Left_End[i]){
              rect(table_out_interval$Left_Start[i]+25,i-0.25,table_out_interval$Left_End[i],i+0.25,border=NA,col="#377eb8BF")
            }
            if (table_out_interval$Right_End[i]-25>=table_out_interval$Right_Start[i]){
              rect(table_out_interval$Right_Start[i],i-0.25,table_out_interval$Right_End[i]-25,i+0.25,border=NA,col="#377eb8BF")
            }
            polygon(x=c(table_out_interval$Right_End[i]-25,table_out_interval$Right_End[i],table_out_interval$Right_End[i]-25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#377eb8BF")
            polygon(x=c(table_out_interval$Left_Start[i]+25,table_out_interval$Left_Start[i],table_out_interval$Left_Start[i]+25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#377eb8BF")
            lines(x=c(table_out_interval$Left_End[i],table_out_interval$Right_Start[i]),y=c(i,i),col="#377eb8BF")
  
          } else if (table_out_interval$Read_Type[i]=="FF"){
            if (table_out_interval$Left_End[i]-25>=table_out_interval$Left_Start[i]){
              rect(table_out_interval$Left_Start[i],i-0.25,table_out_interval$Left_End[i]-25,i+0.25,border=NA,col="#7bccc4BF")
            }
            if (table_out_interval$Right_End[i]-25>=table_out_interval$Right_Start[i]){
              rect(table_out_interval$Right_Start[i],i-0.25,table_out_interval$Right_End[i]-25,i+0.25,border=NA,col="#7bccc4BF")
            }
            polygon(x=c(table_out_interval$Left_End[i]-25,table_out_interval$Left_End[i],table_out_interval$Left_End[i]-25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#7bccc4BF")
            polygon(x=c(table_out_interval$Right_End[i]-25,table_out_interval$Right_End[i],table_out_interval$Right_End[i]-25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#7bccc4BF")
            lines(x=c(table_out_interval$Left_End[i],table_out_interval$Right_Start[i]),y=c(i,i),col="#7bccc4BF")
  
          } else if (table_out_interval$Read_Type[i]=="RR"){ #RR
            if (table_out_interval$Left_Start[i]+25<=table_out_interval$Left_End[i]){
              rect(table_out_interval$Left_Start[i]+25,i-0.25,table_out_interval$Left_End[i],i+0.25,border=NA,col="#7bccc4BF")
            }
            if (table_out_interval$Right_Start[i]+25<=table_out_interval$Right_End[i]){
              rect(table_out_interval$Right_Start[i]+25,i-0.25,table_out_interval$Right_End[i],i+0.25,border=NA,col="#7bccc4BF")
            }
            polygon(x=c(table_out_interval$Left_Start[i]+25,table_out_interval$Left_Start[i],table_out_interval$Left_Start[i]+25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#7bccc4BF")
            polygon(x=c(table_out_interval$Right_Start[i]+25,table_out_interval$Right_Start[i],table_out_interval$Right_Start[i]+25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#7bccc4BF")
            lines(x=c(table_out_interval$Left_End[i],table_out_interval$Right_Start[i]),y=c(i,i),col="#7bccc4BF")
  
          } else if (table_out_interval$Read_Type[i]=="Split-read"){
            rect(table_out_interval$Left_Start[i],i-0.25,table_out_interval$Left_End[i],i+0.25,border=NA,col="#984ea3")
            rect(table_out_interval$Right_Start[i],i-0.25,table_out_interval$Right_End[i],i+0.25,border=NA,col="#984ea3")
            lines(x=c(table_out_interval$Left_End[i],table_out_interval$Right_Start[i]),y=c(i,i),col="#984ea3")
          }
        }
        
        rect(master_info$start,-100,master_info$end,length(table_out_interval$Chr)*2.1,col="#8dd3c71A",border="#8dd3c74D")
        legend("bottomleft",legend=c("-->  <-- [删除]","<--  --> [串联dup或易位]","-->  --> or <--  <-- [反演]","拆分-读取"),
               fill=c("#e41a1c","#377eb8","#7bccc4","#984ea3"),border=NA,cex=0.8,bg="white")
      } else{
        plot(0,0,col=rgb(0,0,0,0),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
        box()
        text(0.5,0.5,"基于过滤器未找到读取。",cex=2)
      }
    }
  })

  observeEvent(input$click_reads,{
    req(rv$config_loaded)
    brush=input$brush_reads
    if (!is.null(brush)){
      ranges$x=c(brush$xmin,brush$xmax)
    }else{
      ranges$x=ranges_default$x
    }
  })
    
  output$reads_table=DT::renderDataTable({
    req(rv$config_loaded)
    if (load$reads > 0) {
      sorted_table <- out_reads$table[sort(out_reads$table$Left_Start,index.return=T,decreasing=F)$ix,]
      # Change the column names for display purposes only (not affecting the data)
      colnames(sorted_table) <- c("读取姓名", "旗标", "Chr", "左起动", "左端", "右启动", "右端", "大小", "MQ", "CIGAR", "读取类型")
      sorted_table
    }
  },options = list(pageLength = 10, language=list(info="显示 _START_ 个条目中的 _END_ 到 _TOTAL_ 个",paginate=list(previous="上一个",`next`="下一个"),search = "搜索:",lengthMenu ="显示 _MENU_ 个条目")),
    rownames = FALSE,
    selection = 'none')
  
  output$gnomadsv_clinvar=renderPlot({
    req(rv$config_loaded, input$cnv_name, master_info$loaded)
    sv1$gnomad=NULL
    sv1$clinvar=NULL
    sv1$cgc=NULL
    
    x9=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script06_file1.txt.gz",sep=""),comment.char="",sep="\t",fill=T)
    x9_gnomad=x9[which(x9$V1=="gnomAD_SV"),]
    
    if (length(as.data.frame(x9_gnomad)[,1])>0){
      x9_gnomad$V4=as.numeric(as.character(x9_gnomad$V4))
      x9_gnomad$V5=as.numeric(as.character(x9_gnomad$V5))
      x9_gnomad$V7=as.numeric(as.character(x9_gnomad$V7))
      gnomad_af_filt=x9_gnomad[c(which(x9_gnomad$V7>=as.numeric(input$gnomad_af)),which(is.na(x9_gnomad$V7)==T)),]
      gnomad_interval=gnomad_af_filt[intersect(which(gnomad_af_filt$V5>=ranges$x[1]),which(gnomad_af_filt$V4<=ranges$x[2])),]
      sv1$gnomad=gnomad_interval
    }
    
    x9_clinvar=x9[which(x9$V1=="ClinVar"),]
    
    if (length(as.data.frame(x9_clinvar)[,1])>0){
      x9_clinvar$V3=as.numeric(as.character(x9_clinvar$V3))
      x9_clinvar$V4=as.numeric(as.character(x9_clinvar$V4))
      clinvar_len=x9_clinvar$V4-x9_clinvar$V3+1
      clinvar_len_filt=x9_clinvar[which(clinvar_len<=as.numeric(input$clinvar_size)),]
      
      if (input$clinvar_subset=="all"){
        clinvar_consequence_filt=clinvar_len_filt
      } else if (input$clinvar_subset=="non-b"){
        clinvar_consequence_filt=clinvar_len_filt[which(clinvar_len_filt$V11!="B"),]
      } else if (input$clinvar_subset=="b"){
        clinvar_consequence_filt=clinvar_len_filt[which(clinvar_len_filt$V11=="B"),]
      } else if (input$clinvar_subset=="p"){
        clinvar_consequence_filt=clinvar_len_filt[which(clinvar_len_filt$V11=="P"),]
      }
      clinvar_interval=clinvar_consequence_filt[intersect(which(clinvar_consequence_filt$V4>=ranges$x[1]),which(clinvar_consequence_filt$V3<=ranges$x[2])),]
      sv1$clinvar=clinvar_interval
    }
    
    x9=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script06_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="",fill=T)
    x9_cgc=x9[which(x9$V1=="cohort_CNV"),]
    
    if (length(as.data.frame(x9_cgc)[,1])>0){
      x9_cgc$V3=as.numeric(as.character(x9_cgc$V3))
      x9_cgc$V4=as.numeric(as.character(x9_cgc$V4))
      cgc_interval=x9_cgc[intersect(which(x9_cgc$V4>=ranges$x[1]),which(x9_cgc$V3<=ranges$x[2])),]
      cgc_interval_sorted=cgc_interval[sort(as.character(cgc_interval[,6]),index.return=T)$ix,]
      sv1$cgc=cgc_interval_sorted
    }
    
    par(mgp=c(2.1,1,0),mai=c(0.4,0.62,0.5,0.1))
    plot(0,0,col=rgb(0,0,0,0),xlim=ranges$x,ylim=c(-20,10),axes=F,xlab="",ylab="",xaxs="i",yaxs="i",main="外部和内部变体数据库")
    box()
    axis(1)
    abline(h=0)
    abline(h=-10)
    
    if(length(sv1$gnomad[,1])>0){
      step=10/(length(sv1$gnomad[,1])+1)
      for(i in 1:length(sv1$gnomad[,1])){
        level=10-step*i
        rect(sv1$gnomad[i,4],level-0.25*step,sv1$gnomad[i,5],level+0.25*step,border="#969696",col="#969696")
        if(length(sv1$gnomad[,1])<=25){
          if (sv1$gnomad[i,5]>=ranges$x[2]){
            text(ranges$x[2],level,paste(sv1$gnomad[i,9],", AF=",sv1$gnomad[i,7],sep=""),cex=0.8,pos=2)
          } else{
            text(sv1$gnomad[i,5],level,paste(sv1$gnomad[i,9],", AF=",sv1$gnomad[i,7],sep=""),cex=0.8,pos=4)
          }
        }
      }
    } else{
      text((ranges$x[1]+ranges$x[2])/2,5,"基于过滤器未发现gnomAD-SV变体。",cex=2)
    }
    
    if (length(sv1$clinvar[,1])>0){
      step=10/(length(sv1$clinvar[,1])+1)
      for(i in 1:length(sv1$clinvar[,1])){
        level=0-step*i
        col_clinvar="#969696"
        if (sv1$clinvar[i,11]=="P"){
          col_clinvar="#f781bf"
        } else if (sv1$clinvar[i,11]=="B"){
          col_clinvar="#e6f5c9"
        }
        
        rect(sv1$clinvar[i,3],level-0.25*step,sv1$clinvar[i,4],level+0.25*step,border=col_clinvar,col=col_clinvar)
        if(length(sv1$clinvar[,1])<=25){
          clinvar_name1=unlist(strsplit(unlist(strsplit(sv1$clinvar[i,10],"<"))[2],">"))[2]
          if (sv1$clinvar[i,4]>=ranges$x[2]){
            text(ranges$x[2],level,clinvar_name1,cex=0.8,pos=2)
          } else{
            text(sv1$clinvar[i,4],level,clinvar_name1,cex=0.8,pos=4)
          }
        }
      }
    } else{
      text((ranges$x[1]+ranges$x[2])/2,-5,"基于过滤器未发现Clinvar变体。",cex=2)
    }
    
    if (length(sv1$cgc[,1])>0){
      step=10/(length(sv1$cgc[,1])+1)
      for(i in 1:length(sv1$cgc[,1])){
        level=-10-step*i
        col_cgc="#969696"
        if (sv1$cgc[i,8]=="1"){
          col_cgc="#fb8072"
        }
        
        rect(sv1$cgc[i,3],level-0.25*step,sv1$cgc[i,4],level+0.25*step,border=col_cgc,col=col_cgc)
        if(length(sv1$cgc[,1])<=25){
          if (sv1$cgc[i,4]>=ranges$x[2]){
            text(ranges$x[2],level,sv1$cgc[i,6],cex=0.7,pos=2)
          } else{
            text(sv1$cgc[i,4],level,sv1$cgc[i,6],cex=0.7,pos=4)
          }
        }
      }
    } else{
      text((ranges$x[1]+ranges$x[2])/2,-15,"基于过滤器未发现其他内部样品中的变体。",cex=2)
    }

    text(ranges$x[1],0.8,"gnomAD-SV",pos=4)
    text(ranges$x[1],-0.8,"ClinVar",pos=4)
    text(ranges$x[1],-10.8,paste("内部队列(合计=",length(sv1$cgc[,1]),")",sep=""),pos=4)
    
    legend("bottomleft",legend=c("ClinVar P/LP","ClinVar B/LB","同一个家庭"),fill=c("#f781bf","#e6f5c9","#fb8072"),cex=0.8,bg="white",border=NA)
    rect(master_info$start,-1000,master_info$end,1000,col="#8dd3c71A",border="#8dd3c74D")
  })

  observeEvent(input$click_gnomadsv_clinvar,{
    req(rv$config_loaded)
    brush=input$brush_gnomadsv_clinvar
    if (!is.null(brush)){
      ranges$x=c(brush$xmin,brush$xmax)
    }else{
      ranges$x=ranges_default$x
    }
  })
  
  output$gnomadsv_table=DT::renderDataTable({
    req(rv$config_loaded)
    if(is.null(sv1$gnomad)==F){
      x10=sv1$gnomad[,c(3:8,2)]
      colnames(x10)=c("Chr","开始","结束","类型","Popmax AF","QC标志","ID")
      x10
    } else{
      datatable(as.data.frame("基于过滤器未发现gnomAD-SV变体。"),selection='none',options=list(dom="t",ordering=F,
                                                                                                       headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F)
    }
  },escape=F,rownames=F,selection='none',options = list(pageLength = 10, language=list(info="显示 _START_ 个条目中的 _END_ 到 _TOTAL_ 个",paginate=list(previous="上一个",`next`="下一个"),search = "搜索:",lengthMenu ="显示 _MENU_ 个条目")))
  
  output$clinvar_table=DT::renderDataTable({
    req(rv$config_loaded)
    if(is.null(sv1$clinvar)==F){
      x11=sv1$clinvar[,c(2:4,6,5,7:10)]
      colnames(x11)=c("Chr","开始","结束","类型","口译","条件","等位基因起源","基因","ID")
      x11
    } else{
      datatable(as.data.frame("基于过滤器未发现Clinvar变体。"),selection='none',options=list(dom="t",ordering=F,
                                                                                                       headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F)
    }
  },escape=F,rownames=F,selection='none',options = list(pageLength = 10, language=list(info="显示 _START_ 个条目中的 _END_ 到 _TOTAL_ 个",paginate=list(previous="上一个",`next`="下一个"),search = "搜索:",lengthMenu ="显示 _MENU_ 个条目")))
  
  output$cgc_table=DT::renderDataTable({
    if(is.null(sv1$cgc)==F){
      x12=sv1$cgc[,c(2:7)]
      colnames(x12)=c("Chr","开始","结束","类型","样品","算法")
      x12
    } else{
      datatable(as.data.frame("基于过滤器未发现其他CGC样品中的变体。"),selection='none',options=list(dom="t",ordering=F,
                                                                                                     headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F)
    }
  },escape=F,rownames=F,selection='none',options = list(pageLength = 10, language=list(info="显示 _START_ 个条目中的 _END_ 到 _TOTAL_ 个",paginate=list(previous="上一个",`next`="下一个"),search = "搜索:",lengthMenu ="显示 _MENU_ 个条目")))
  
  output$genes_table=DT::renderDataTable({
    req(rv$config_loaded, input$cnv_name, master_info$loaded)
    x13=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script07_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
    x13=x13[,c(1,19,5:17)]
    colnames(x13)=c("姓名","钢绞线","重叠","HI","TS","pLI","LOEUF","pRec","pNull","表达式","GO","OMIM","GenCC","链接","搜索")
    
    datatable(x13,escape=F,class="cell-border stripe",rownames=NULL,options=list(pageLength=25,language=list(info="显示 _START_ 个条目中的 _END_ 到 _TOTAL_ 个",paginate=list(previous="上一个",`next`="下一个"),search = "搜索:",lengthMenu ="显示 _MENU_ 个条目")),selection=list(mode="single", target="cell")) %>%
      DT::formatStyle(columns=12:13,fontSize='90%') %>%
      DT::formatStyle(columns=c(6:10,14,15),"text-align"='center') %>%
      DT::formatStyle('pLI',backgroundColor=styleInterval(c(0.5,0.9),c("white","#ffffb2","#fcae91"))) %>%
      DT::formatStyle('LOEUF',backgroundColor=styleInterval(c(0.35),c("#fcae91","white"))) %>%
      DT::formatStyle('pRec',backgroundColor=styleInterval(c(0.9),c("white","#bdc9e1"))) %>%
      DT::formatStyle(columns=c(15),width='100px')
  })

  observeEvent(input$genes_table_cells_selected,{
    req(rv$config_loaded, input$cnv_name, master_info$loaded)
    x13_v1=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script07_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
    coord=unlist(input$genes_table_cells_selected)
    if (length(coord)==2){
      if (coord[2]==1-1 || coord[2]==2-1 || coord[2]==3-1){
        showModal(mymodal(title="成绩单信息",easyClose=T,
                          HTML(x13_v1$V20[as.numeric(coord[1])]),idcss="xl1"))
      }
    }
  })
    
  output$dosage_table=DT::renderDataTable({
    req(rv$config_loaded, input$cnv_name, master_info$loaded)
    x9=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script06_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="",fill=T)
    x9_dosage=x9[which(x9$V1=="ClinGen_Region"),]
    hi_region_index=c(which(x9_dosage$V6=="Sufficient (3)"),which(x9_dosage$V6=="Emerging (2)"),which(x9_dosage$V6=="Little (1)"),
                      which(x9_dosage$V6=="Recessive (30)"),which(x9_dosage$V6=="Unlikely (40)"))
    ts_region_index=c(which(x9_dosage$V7=="Sufficient (3)"),which(x9_dosage$V7=="Emerging (2)"),which(x9_dosage$V7=="Little (1)"),
                      which(x9_dosage$V7=="Recessive (30)"),which(x9_dosage$V7=="Unlikely (40)"))

    ds_region_index=union(hi_region_index,ts_region_index)
    if (master_info$type=="DEL"){
      ds_region_index=hi_region_index
    }
    if (master_info$type=="DEL"){
      ds_table=x9_dosage[ds_region_index,2:6]
      colnames(ds_table)=c("姓名","Chr","开始","结束","Haploinsufficiency")
    } else{
      ds_table=x9_dosage[ds_region_index,2:7]
      colnames(ds_table)=c("姓名","Chr","开始","结束","Haploinsufficiency","Triplosensitivity")      
    }
    ds_table
  },escape=F,rownames=F,options=list(pageLength=25,language=list(info="显示 _START_ 个条目中的 _END_ 到 _TOTAL_ 个",paginate=list(previous="上一个",`next`="下一个"),search = "搜索:",lengthMenu ="显示 _MENU_ 个条目")),selection='none')
  
  output$genes=renderPlot({
    req(rv$config_loaded, input$cnv_name, master_info$loaded)
    par(mgp=c(2.1,1,0),mai=c(0.4,0.62,0.5,0.1))
    
    plot(0,0,col=rgb(0,0,0,0),
         xlab="",xlim=ranges$x,ylab="",ylim=c(0,3.5),xaxs="i",yaxs="i",axes=F,main="基因组邻域")
    box()
    
    if (load$genes==0){
      text((as.numeric(ranges$x[1])+as.numeric(ranges$x[2]))/2,1.75,"单击上面的 \"绘制基因组邻域\" 按钮以加载此图。",cex=2)
    }
    if (load$genes>0){
      x9=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script06_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="",fill=T)
      x9_dosage=x9[which(x9$V1=="ClinGen_Region"),]
      if (length(as.data.frame(x9_dosage)[,1])>0){
        hi_region_index=c(which(x9_dosage$V6=="Sufficient (3)"),which(x9_dosage$V6=="Emerging (2)"),which(x9_dosage$V6=="Little (1)"),
                          which(x9_dosage$V6=="Recessive (30)"),which(x9_dosage$V6=="Unlikely (40)"))
        ts_region_index=c(which(x9_dosage$V7=="Sufficient (3)"),which(x9_dosage$V7=="Emerging (2)"),which(x9_dosage$V7=="Little (1)"),
                          which(x9_dosage$V7=="Recessive (30)"),which(x9_dosage$V7=="Unlikely (40)"))
        ds_region_index=union(hi_region_index,ts_region_index)
        if (master_info$type=="DEL"){
          ds_region_index=hi_region_index
        }
        sv1$dosage_hi=x9_dosage[hi_region_index,]
        sv1$dosage_ts=x9_dosage[ts_region_index,]
      }
  
      axis(1)
      axis(2,at=seq(0,1,0.2))
      abline(h=c(1,2,2.5,3))
      
      if (length(sv1$dosage_hi[,1])>0){
        for (i in 1:length(sv1$dosage_hi[,1])){
          col_hi="#969696"
          if (sv1$dosage_hi[i,6]=="Sufficient (3)"){
            col_hi="#f03b20"
          } else if (sv1$dosage_hi[i,6]=="Emerging (2)"){
            col_hi="#fd8d3c"
          } else if (sv1$dosage_hi[i,6]=="Little (1)"){
            col_hi="#fecc5c"
          }
          
          step_hi=0.5/(length(sv1$dosage_hi[,1])+1)
          rect(sv1$dosage_hi[i,4],2.5+(i-0.25)*step_hi,sv1$dosage_hi[i,5],2.5+(i+0.25)*step_hi,border=NA,col=col_hi)
          hi_name1=unlist(strsplit(unlist(strsplit(sv1$dosage_hi[i,2],"<"))[2],">"))[2]
          if (sv1$dosage_hi[i,5]>=ranges$x[2]){
            if(sv1$dosage_hi[i,4]<=ranges$x[2]){      
              text(ranges$x[2],2.5+i*step_hi,hi_name1,cex=0.7,pos=2)
            }
          } else{
            text(sv1$dosage_hi[i,5],2.5+i*step_hi,hi_name1,cex=0.7,pos=4)
          }
        }
      } else{
        text((ranges$x[1]+ranges$x[2])/2,2.75,"此变体与Clingen HI区域不重叠。")
      }
      
      if (master_info$type=="DEL"){
        text((ranges$x[1]+ranges$x[2])/2,3.25,"TS不适用，因为变体是DEL。")
      } else{
        if (length(sv1$dosage_ts[,1])>0){
          for (i in 1:length(sv1$dosage_ts[,1])){
            col_ts="#969696"
            if (sv1$dosage_ts[i,7]=="Sufficient (3)"){
              col_ts="#f03b20"
            } else if (sv1$dosage_ts[i,7]=="Emerging (2)"){
              col_ts="#fd8d3c"
            } else if (sv1$dosage_ts[i,7]=="Little (1)"){
              col_ts="#fecc5c"
            }
          
            step_ts=0.5/(length(sv1$dosage_ts[,1])+1)
            rect(sv1$dosage_ts[i,4],3+(i-0.25)*step_ts,sv1$dosage_ts[i,5],3+(i+0.25)*step_ts,border=NA,col=col_ts)
            ts_name1=unlist(strsplit(unlist(strsplit(sv1$dosage_ts[i,2],"<"))[2],">"))[2]
            if (sv1$dosage_ts[i,5]>=ranges$x[2]){
              if(sv1$dosage_ts[i,4]<=ranges$x[2]){ 
                text(ranges$x[2],3+i*step_ts,ts_name1,cex=0.7,pos=2)
              }
            } else{
              text(sv1$dosage_ts[i,5],3+i*step_ts,ts_name1,cex=0.7,pos=4)
            }
          }
        } else{
          text((ranges$x[1]+ranges$x[2])/2,3.25,"此变体不与Clingen TS区域重叠。")
        }
      }
      
      x15=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script05_file3.txt.gz",sep=""),comment.char="",sep="\t",quote="")
      y_step=1/(max(x15$V2)+2)
      for (i in 1:length(x15$V1)){
        level=1+y_step*(x15$V2[i]+1)
        pos=unlist(strsplit(x15[i,3],"\\:"))
        for (j in seq(1,length(pos),2)){
          rect(pos[j],level-y_step*0.25,pos[j+1],level+y_step*0.25,border="black",col="black")
          if (j!=1){
            lines(x=c(pos[j-1],pos[j]),y=c(level,level))
          }
          if (ranges$x[2]-ranges$x[1]<=5e4){
            if ((j+1)%%4==0){
              text((as.numeric(pos[j])+as.numeric(pos[j+1]))/2,level+y_step*0.25,x15$V1[i],pos=3,cex=0.7)
            } else{
              text((as.numeric(pos[j])+as.numeric(pos[j+1]))/2,level-y_step*0.25,x15$V1[i],pos=1,cex=0.7)
            }
          }
        }
      }
      
      x14=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script05_file1.txt.gz",sep=""))
      x14sub1=x14$V1[is.na(x14$V4)==F]
      x14sub4=x14$V4[is.na(x14$V4)==F]
      
      # speed up using plotH
      # for(i in 1:length(x14sub1)){
      #   if (x14sub4[i]>=0.2){
      #     rect(x14sub1[i]-0.5,0,x14sub1[i]+0.5,x14sub4[i],border="#756bb1",col="#756bb1")
      #   } else{
      #     rect(x14sub1[i]-0.5,0,x14sub1[i]+0.5,x14sub4[i],border="#969696",col="#969696")        
      #   }
      # }
      par(new=T)
      plotH(x14sub1,x14sub4,width=1,col="#756bb1",xlab="",ylab="",xlim=ranges$x,ylim=c(0,3.5),xaxs="i",yaxs="i",axes=F,border="#756bb1")

      x17=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script07_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
      gnomad_pli_gene=c()
      gnomad_pli_start=c()
      gnomad_pli_end=c()
      gnomad_pli_cat=c()
      for (i in 1:length(x17$V1)){
        if(is.na(x17$V8[i])==F){
          if (x17$V8[i]>=0.9 || x17$V9[i]<=0.35){
            gnomad_pli_gene=c(gnomad_pli_gene,unlist(strsplit(unlist(strsplit(x17$V1[i],"<"))[2],">"))[2])
            gnomad_pli_start=c(gnomad_pli_start,x17$V3[i])
            gnomad_pli_end=c(gnomad_pli_end,x17$V4[i])
            gnomad_pli_cat=c(gnomad_pli_cat,1)
          } else if (x17$V8[i]>=0.5){
            gnomad_pli_gene=c(gnomad_pli_gene,unlist(strsplit(unlist(strsplit(x17$V1[i],"<"))[2],">"))[2])
            gnomad_pli_start=c(gnomad_pli_start,x17$V3[i])
            gnomad_pli_end=c(gnomad_pli_end,x17$V4[i])
            gnomad_pli_cat=c(gnomad_pli_cat,2)
          } else if (x17$V8[i]>=0.1){
            gnomad_pli_gene=c(gnomad_pli_gene,unlist(strsplit(unlist(strsplit(x17$V1[i],"<"))[2],">"))[2])
            gnomad_pli_start=c(gnomad_pli_start,x17$V3[i])
            gnomad_pli_end=c(gnomad_pli_end,x17$V4[i])
            gnomad_pli_cat=c(gnomad_pli_cat,3)
          }
        }
      }

      if (length(gnomad_pli_cat)>0){
        gnomad_pli_color=rep("#ffffcc",length(gnomad_pli_cat))
        gnomad_pli_color[which(gnomad_pli_cat==2)]="#a1dab4"
        gnomad_pli_color[which(gnomad_pli_cat==1)]="#2c7fb8"
        
        x17_step=0.5/(length(gnomad_pli_cat)+1)
        for (i in 1:length(gnomad_pli_cat)){
          rect(gnomad_pli_start[i],2+(i-0.25)*x17_step,gnomad_pli_end[i],2+(i+0.25)*x17_step,border=gnomad_pli_color[i],col=gnomad_pli_color[i])
          if (gnomad_pli_end[i]>=ranges$x[2]){
            if(gnomad_pli_start[i]<=ranges$x[2]){      
              text(ranges$x[2],2+i*x17_step,gnomad_pli_gene[i],cex=0.7,pos=2)
            }
          } else{
            text(gnomad_pli_end[i],2+i*x17_step,gnomad_pli_gene[i],cex=0.7,pos=4)
          }
        }
      } else{
        text((as.numeric(ranges$x[1])+as.numeric(ranges$x[2]))/2,2.25,"该变体与pLI >= 0.1的基因不重叠。")
      }
            
      text(ranges$x[1],3.4,"ClinGen Triplosensitivity Map",pos=4)
      text(ranges$x[1],2.9,"ClinGen Haploinsufficiency Map",pos=4)
      text(ranges$x[1],2.4,"gnomAD HI 约束",pos=4)
      text(ranges$x[1],1.9,"基因",pos=4)
      text(ranges$x[1],0.9,"pext",pos=4)
      legend("topright",legend=c("足够的","新兴","小","隐性/不太可能","pLI 0.9+ / LOEUF 0.35-","pLI 0.5+","pLI 0.1+"),
             fill=c("#f03b20","#fd8d3c","#fecc5c","#969696","#2c7fb8","#a1dab4","#ffffcc"),cex=0.8,bg="white",border=NA)
      # legend("bottomleft",legend=c("pext >= 0.2"),fill=c("#756bb1"),cex=0.8,bg="white",border=NA) # currently all pext bars (regardless of level) are in purple
      rect(master_info$start,-1000,master_info$end,1000,col="#8dd3c71A",border="#8dd3c74D")
    }
  })
  
  observeEvent(input$click_genes,{
    req(rv$config_loaded)
    brush=input$brush_genes
    if (!is.null(brush)){
      ranges$x=c(brush$xmin,brush$xmax)
    }else{
      ranges$x=ranges_default$x
    }
  })
  
  output$version_control1=renderTable({
    req(rv$config_loaded, input$cnv_name, master_info$loaded)
    x16=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".log.txt.gz",sep=""),fill=T,comment.char="",sep="\t",quote="",head=F)
    x16[1:3,1:2]
  },colnames=F,spacing="s",striped=T,bordered=T)

  output$version_control2=renderTable({
    req(rv$config_loaded, input$cnv_name, master_info$loaded)
    x16=read.table(paste(rv$temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".log.txt.gz",sep=""),fill=T,comment.char="",sep="\t",quote="",head=F)
    x16[4:length(x16$V1),1:3]
  },colnames=F,spacing="s",striped=T,bordered=T)
}

shinyApp(ui=ui,server=server,options = list(port = 4704, host = "0.0.0.0"))
