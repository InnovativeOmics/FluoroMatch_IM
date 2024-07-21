##################################################
## PROJECT: Ion Mobility
## SCRIPT:  Scoring
## DATE:    2024-07-14
## AUTHORs:  Paul Stelben, Jeremy Koelmel, Michael Kummer, David Schiessel
##################################################

# Clear environment
rm(list = ls())
closeAllConnections()

ManuallyInputVariables <- FALSE
RT_flagging <- TRUE #JPK: for PFAS analysis
ParallelComputing <- TRUE
TWeen_pos <- FALSE #PJS: for PolyMatch

if (ParallelComputing == TRUE) {
  if("foreach" %in% rownames(installed.packages()) == FALSE) {install.packages("foreach")}
  library(foreach)
  if("doParallel" %in% rownames(installed.packages()) == FALSE) {install.packages("doParallel")}
  library(doParallel)
  if("parallelly" %in% rownames(installed.packages()) == FALSE) {install.packages("parallelly", repos = "http://cran.us.r-project.org")}
  library("parallelly")
  ###NOTE not sure if any of this is necessary but was causing error
  DC <- as.numeric(availableCores(constraints = "connections"))
  print(paste0(freeConnections()," connections are free and will be used of ",availableConnections()," R hard limited connections; the following numbers of cores exist on your system: ",as.numeric(detectCores())))
  if (is.na(DC)) {
    DC <- makePSOCKcluster(4)
    registerDoParallel(DC)
  } else if (DC <= 4) {
    DC <- makePSOCKcluster(DC)
    registerDoParallel(DC)
  } else {
    DC <- makePSOCKcluster(DC-2)
    registerDoParallel(DC)
  }
  # Force the directory for parallel computing to be the one that is distributed not the local R directory (doPar error)
  invisible(clusterEvalQ(DC, .libPaths()[2]))
}

# # if("tictoc" %in% rownames(installed.packages()) == FALSE) {install.packages("tictoc")}
# # library(tictoc)
# # tic.clear()
# # tic("Main")
# if("sqldf" %in% rownames(installed.packages()) == FALSE) {install.packages("sqldf")}
if("Rdisop" %in% rownames(installed.packages()) == FALSE) {install.packages("Rdisop")}
# if("RSQLite" %in% rownames(installed.packages()) == FALSE) {install.packages("RSQLite")}
if("gWidgets" %in% rownames(installed.packages()) == FALSE) {install.packages("gWidgets")}
if("gWidgetstcltk" %in% rownames(installed.packages()) == FALSE) {install.packages("gWidgetstcltk")}
require(gWidgets)
require(gWidgetstcltk)
options(guiToolkit="tcltk")
if("tcltk2" %in% rownames(installed.packages()) == FALSE) {install.packages("tcltk2")}
library(Rdisop)
# library(RSQLite)
# library(sqldf)
options(warn=-1)#suppress warning on
if (!require("stringr", quietly = TRUE))install.packages("stringr", repos = "http://cran.us.r-project.org")
library("stringr")

##################################################
## SECTION 1: Read in CSVs and input variables
##################################################

if (ManuallyInputVariables == TRUE) {
  # IDedTable_dir <- "C:/Users/pstel/Documents/PFAS/2023_04_19_IonMobility/Input_Hannah_Florance_Surfactants_Neg.csv"
  IDedTable_dir <- "C:/Users/pstel/Documents/PFAS/2023_04_19_IonMobility/Example_Data/CM_PFAS_Level6_IM4_bit_DI3_d_DeMP_Neg.csv"
  Library_dir <- "C:/Users/pstel/Documents/PFAS/2023_04_19_IonMobility/2023_Weil_PFAS_CCS_Libs.csv"
  RepeatingUnits_dir <- "C:/Users/pstel/Documents/PFAS/2023_04_19_IonMobility/REPEATING_UNITS_INPUT.csv"
  
  # IDedTable columns
  Mass_col <- 5
  Retention_col <- 2
  CCS_col <- 4
  DT_col <- 3
  rowID_col <- 1
  
  # Library columns - label as 0 if not present
  Name_col <- 1
  ID_col <- 0
  Form_col <- 2
  SMILES_col <- 0
  
  RowStartForFeatureTableData <- 1
  RT_flagging <- TRUE
  DT_flagging <- TRUE
  mz_win <- .01 # PrecursorMassAccuracy
  mz_ppm_win <- 10  # Precursor Mass Accuracy in  ppm
  CCS_per <- 2
  upper <- 0.12
  lower <- -0.11
  NegPos <- "Neg"
  
} else {
  done <- tclVar(0)
  
  Mass_col <- NULL
  Retention_col <- NULL
  CCS_col <- NULL
  DT_col <- NULL
  rowID_col <- NULL
  mz_win <- NULL
  mz_ppm_win <- NULL
  Name_col <- NULL
  ID_col <- NULL
  Form_col <- NULL
  SMILES_col <- NULL
  
  RowStartForFeatureTableData <- NULL
  RT_flagging <- FALSE
  DT_flagging <- FALSE

  CCS_per <- NULL
  upper <- NULL
  lower <- NULL
  NegPos <- NULL
  
  GUILauncher <- function() {
    
    Mass_col <<- tclVar("5")
    Retention_col <<- tclVar("2")
    CCS_col <<- tclVar("4")
    DT_col <<- tclVar("3")
    rowID_col <<- tclVar("1")
    
    Name_col <<- tclVar("1")
    ID_col <<- tclVar("0")
    Form_col <<- tclVar("2")
    SMILES_col <<- tclVar("23")
    
    RowStartForFeatureTableData <<- tclVar("2")
    mz_win <<- tclVar(".01")
    mz_ppm_win <<- tclVar("10")
    CCS_per <<-tclVar("2")
    upper <<- tclVar("0.12")
    lower <<- tclVar("-0.11")
    NegPos <<- tclVar("Neg")
    
    tt <- tktoplevel()
    tkwm.title(tt,"MAIN_ALL Settings")
    
    Mass_col.entry <- tkentry(tt, textvariable = Mass_col)
    Retention_col.entry <- tkentry(tt, textvariable = Retention_col)
    CCS_col.entry <- tkentry(tt, textvariable = CCS_col)
    DT_col.entry <- tkentry(tt, textvariable = DT_col)
    rowID_col.entry <- tkentry(tt, textvariable = rowID_col)
    
    Name_col.entry <- tkentry(tt, textvariable = Name_col)
    ID_col.entry <- tkentry(tt, textvariable = ID_col)
    Form_col.entry <- tkentry(tt, textvariable = Form_col)
    SMILES_col.entry <- tkentry(tt, textvariable = SMILES_col)
    
    RowStartForFeatureTableData.entry <- tkentry(tt, textvariable = RowStartForFeatureTableData)
    mz_win.entry <- tkentry(tt, textvariable = mz_win)
    mz_ppm_win.entry <- tkentry(tt, textvariable = mz_ppm_win)
    CCS_per.entry <- tkentry(tt, textvariable = CCS_per)
    upper.entry <- tkentry(tt, textvariable = upper)
    lower.entry <- tkentry(tt, textvariable = lower)
    NegPos.entry <- tkentry(tt, textvariable = NegPos)
    
    submit <- function() {
      Mass_col <<- as.numeric(tclvalue(Mass_col))
      Retention_col <<- as.numeric(tclvalue(Retention_col))
      CCS_col <<- as.numeric(tclvalue(CCS_col))
      DT_col <<- as.numeric(tclvalue(DT_col))
      rowID_col <<- as.numeric(tclvalue(rowID_col))
      
      Name_col <<- as.numeric(tclvalue(Name_col))
      ID_col <<- as.numeric(tclvalue(ID_col))
      Form_col <<- as.numeric(tclvalue(Form_col))
      SMILES_col <<- as.numeric(tclvalue(SMILES_col))
      
      RowStartForFeatureTableData <<- as.numeric(tclvalue(RowStartForFeatureTableData))
      mz_win <<- as.numeric(tclvalue(mz_win))
      mz_ppm_win <<- as.numeric(tclvalue(mz_ppm_win))
      CCS_per <<- as.numeric(tclvalue(CCS_per))
      upper <<- as.numeric(tclvalue(upper))
      lower <<- as.numeric(tclvalue(lower))
      NegPos <<- as.character(tclvalue(NegPos))
      
      tclvalue(done) <- 1
      tkdestroy(tt)
    }
    
    submit.but <- tkbutton(tt, text="Run", command=submit)
    
    
    IDedTable_fun <- function() {
      IDedTable_dir <<- tk_choose.files(default = "", caption = "IDedTable file", multi = FALSE, filters = NULL, index = 1)
    }
    IDedTable_fun.but <- tkbutton(tt, text="Select File", command = IDedTable_fun)
    
    Library_fun <- function() {
      Library_dir <<- tk_choose.files(default = "", caption = "Library file", multi = FALSE, filters = NULL, index = 1)
    }
    Library_fun.but <- tkbutton(tt, text="Select File", command = Library_fun)
    
    RepeatingUnits_fun <- function() {
      RepeatingUnits_dir <<- tk_choose.files(default = "", caption = "RepeatingUnits file", multi = FALSE, filters = NULL, index = 1)
    }
    RepeatingUnits_fun.but <- tkbutton(tt, text="Select File", command = RepeatingUnits_fun)
    
    RT_fun <- function() {
      RT_flagging <<- TRUE
    }
    RT_fun.but <- tkbutton(tt, text="Yes", command = RT_fun)
    
    RT_fun_neg <- function() {
      RT_flagging <<- FALSE
    }
    RT_fun_neg.but <- tkbutton(tt, text="No", command = RT_fun_neg)
    
    DT_fun <- function() {
      DT_flagging <<- TRUE
    }
    DT_fun.but <- tkbutton(tt, text="Yes", command = DT_fun)
    
    DT_fun_neg <- function() {
      DT_flagging <<- FALSE
    }
    DT_fun_neg.but <- tkbutton(tt, text="No", command = DT_fun_neg)
    
    tkgrid(tklabel(tt, text="Choose input file with button") , IDedTable_fun.but)
    tkgrid(tklabel(tt, text="Choose library file with button") , Library_fun.but)
    tkgrid(tklabel(tt, text="Choose repeating units file with button") , RepeatingUnits_fun.but)
    tkgrid(tklabel(tt, text=""))
    tkgrid(tklabel(tt, text="Feature Table columns"))
    tkgrid(tklabel(tt,text="Mass_col"), Mass_col.entry)
    tkgrid(tklabel(tt,text="Retention_col"), Retention_col.entry)
    tkgrid(tklabel(tt,text="CCS_col"), CCS_col.entry)
    tkgrid(tklabel(tt,text="DT_col"), DT_col.entry)
    tkgrid(tklabel(tt,text="rowID_col"), rowID_col.entry)
    tkgrid(tklabel(tt, text=""))
    tkgrid(tklabel(tt, text="CCS Library columns"))
    tkgrid(tklabel(tt,text="Name of Analytes (Col)"), Name_col.entry)
    tkgrid(tklabel(tt,text="Identifier (Numeric) (can leave as 0)"), ID_col.entry)
    tkgrid(tklabel(tt,text="Neutral Molecular Formula (Col)"), Form_col.entry)
    tkgrid(tklabel(tt,text="SMILES_col"), SMILES_col.entry)
    tkgrid(tklabel(tt, text=""))
    tkgrid(tklabel(tt, text="Other"))
    tkgrid(tklabel(tt,text="Row Starts For Feature Table Data"), RowStartForFeatureTableData.entry)
    tkgrid(tklabel(tt, text="RT Order flagging (Series)") , RT_fun.but, RT_fun_neg.but)
    tkgrid(tklabel(tt, text="DT Order flagging (Series") , DT_fun.but, DT_fun_neg.but)
    tkgrid(tklabel(tt,text="Mass Window (DA)"), mz_win.entry)
    tkgrid(tklabel(tt,text="Mass Window (ppm)"), mz_ppm_win.entry)
    tkgrid(tklabel(tt,text="CCS Window (%)"), CCS_per.entry)
    tkgrid(tklabel(tt,text="mass defect flag: upper"), upper.entry)
    tkgrid(tklabel(tt,text="mass defect flag: lower"), lower.entry)
    tkgrid(tklabel(tt,text="Polarity ('Neg' or 'Pos')"), NegPos.entry)
    tkgrid(tklabel(tt, text=""))
    tkgrid(submit.but, pady = 8, columnspan=2)
    
  }
  
  GUILauncher()
  tkwait.variable(done)
}

#Import tables
IDedTable <- as.matrix(read.csv(IDedTable_dir, sep = ","))
Library <- as.matrix(read.csv(Library_dir, sep = ","))

if (nrow(IDedTable) < 2) {
  message("No IonMobility data, stopping scoring function")
  stop()
}

# Change name of columns for PowerBI consistency
colnames(IDedTable)[CCS_col] <- "CCS"
colnames(IDedTable)[DT_col] <- "DT"
colnames(IDedTable)[rowID_col] <- "row.ID"

##################################################
## SECTION 2: Identify chemicals based on library
##################################################

ID_Ranked <- matrix("5_NoID", nrow(IDedTable), 1)
colnames(ID_Ranked) <- "ID_Ranked"
#JPK Technically the data starts in row 2 but it looks like you uploaded data with Column Headers on this time changing things
RowStartForFeatureTableData<-RowStartForFeatureTableData-1
PMA <- mz_win / 2
# CCS_win <- PMA # Should be higher??
for (i in 1:nrow(Library)) {
  lib_mz <- as.numeric(Library[i, 4])
  lib_CCS <- as.numeric(Library[i, 14])
  mz_list <- as.numeric(IDedTable[RowStartForFeatureTableData:nrow(IDedTable), Mass_col])
  CCS_list <- as.numeric(IDedTable[RowStartForFeatureTableData:nrow(IDedTable), CCS_col])
  mz_matches <- which(mz_list >= (lib_mz - PMA) & mz_list <= (lib_mz + PMA))
  if (length(mz_matches) > 0) {
    # message("DOG")
    for (j in 1:length(mz_matches)) {
      x <- mz_matches[j]
      if (length(grep("5_", ID_Ranked[x,])) > 0) {
        ID_Ranked[x,] <- paste(paste0("4_", Library[i, Name_col]), Library[i, ID_col], Library[i, Form_col], Library[i, SMILES_col], sep = ";")
      } else if (length(grep("4_", ID_Ranked[x,])) > 0) {
        ID_Ranked[x,] <- paste(ID_Ranked[x,], "|", paste(paste0("4_", Library[i, Name_col]), Library[i, ID_col], Library[i, Form_col], Library[i, SMILES_col], sep = ";"))
      } else if (length(grep("1_", ID_Ranked[x,])) > 0) {
        next
      }
    }
    CCS_win <- ((CCS_per / 100) * lib_CCS) / 2
    CCS_matches <- which(CCS_list[mz_matches] >= (lib_CCS - CCS_win) & CCS_list[mz_matches] <= (lib_CCS + CCS_win))
    if (length(CCS_matches > 0)) {
      # message("CAT")
      for (j in 1:length(CCS_matches)) {
        x <- mz_matches[CCS_matches[j]]
        if (length(grep("5_", ID_Ranked[x,])) > 0) {
          ID_Ranked[x,] <- paste(paste0("1_", Library[i, Name_col]), Library[i, ID_col], Library[i,Form_col], Library[i, SMILES_col], sep = ";")
        } else if (length(grep("4_", ID_Ranked[x,])) > 0) {
          ID_Ranked[x,] <- paste(paste0("1_", Library[i, Name_col]), Library[i, ID_col], Library[i,Form_col], Library[i, SMILES_col], sep = ";")
        } else if (length(grep("1_", ID_Ranked[x,])) > 0) {
          ID_Ranked[x,] <- paste(ID_Ranked[x,], "|", paste(paste0("1_", Library[i, Name_col]), Library[i, ID_col], Library[i,Form_col], Library[i, SMILES_col], sep = ";"))
        }
      }
    }
  }
}

ID_Ranked <- gsub(";;", ";", ID_Ranked)
ID_Ranked <- gsub("; \\|", " \\|", ID_Ranked)
ID_Ranked <- gsub(";$", "", ID_Ranked)

IDedTable <- cbind(IDedTable, ID_Ranked)

##################################################
## SECTION 3: Kendrick Mass Defect
##################################################

#Import RepeatingUnits
RepeatingUnits <- read.csv(RepeatingUnits_dir,sep=",")
RepeatingUnits <- as.matrix(RepeatingUnits)
Names <- rev(RepeatingUnits[which(RepeatingUnits[,3] == TRUE | RepeatingUnits[,3] == "TRUE" | RepeatingUnits[,3] == " TRUE"),1])
Series <- rev(RepeatingUnits[which(RepeatingUnits[,3] == TRUE | RepeatingUnits[,3] == "TRUE" | RepeatingUnits[,3] == " TRUE"),2])

# FOR DEBUGGING
# IDedTable <- IDedTable[1:3000,]

ID_Ranked_col <- which(colnames(IDedTable) == "ID_Ranked")

# add columns and calculate mass defect and Kendrick mass defect for Mass_col
m <- matrix(0, nrow(IDedTable), 2)
IDedTable <- cbind(IDedTable, m)
nrow_IDedTable <- nrow(IDedTable)
ncol_IDedTable <- ncol(IDedTable)
colnames(IDedTable)[ncol_IDedTable - 1] <- "nominal mass"
colnames(IDedTable)[ncol_IDedTable] <- "mass defect"
Mass_vec <- as.numeric(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, Mass_col])
nmass <- round(Mass_vec)
IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 1] <- as.numeric(nmass)
massd <- Mass_vec - nmass
IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable] <- as.numeric(massd)
md_col <- ncol_IDedTable
hseries_count_cols <- c()

# timestamp()

length_Series <- length(Series)
for(x in 1:length_Series){
  # print(x)
  # x <- 1
  Mass_vec <- as.numeric(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, Mass_col])
  exactMass <- Rdisop::getMolecule(Series[x], maxisotopes = 1)$exactmass ##added Rdisop:: as Rcdk has the same function
  Mass <- round(exactMass)
  m <- matrix(0, nrow(IDedTable), 5)
  IDedTable <- cbind(IDedTable, m)
  nrow_IDedTable <- nrow(IDedTable)
  ncol_IDedTable <- ncol(IDedTable)
  colnames(IDedTable)[ncol_IDedTable - 4] <- paste(Names[x], "Kendrick_mass", sep="_")
  colnames(IDedTable)[ncol_IDedTable - 3] <- paste(Names[x], "nominal_Kendrick_mass", sep="_")
  colnames(IDedTable)[ncol_IDedTable - 2] <- paste(Names[x], "Kendrick_mass_defect", sep="_")
  colnames(IDedTable)[ncol_IDedTable - 1] <- paste(Names[x], "KMD_Group", sep="_")
  colnames(IDedTable)[ncol_IDedTable] <- paste(Names[x], "homologous_series", sep="_")
  
  KM <- Mass_vec * (Mass/exactMass)
  IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 4] <- as.numeric(KM)
  nKM <- round(KM)
  IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 3] <- as.numeric(nKM)
  KMD <- KM - nKM
  IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 2] <- as.numeric(KMD)
  
  # Sort and group by KMD
  tmptable <- IDedTable[RowStartForFeatureTableData:nrow_IDedTable,]
  tmptable <- tmptable[order(as.numeric(tmptable[, ncol_IDedTable - 2])),]
  IDedTable[RowStartForFeatureTableData:nrow_IDedTable,] <- tmptable
  KMD_col <- ncol_IDedTable - 2
  group <- 1
  start <- RowStartForFeatureTableData
  IDedTable[RowStartForFeatureTableData, KMD_col + 1] <- group
  for (i in (RowStartForFeatureTableData + 1):nrow_IDedTable) {
    if (as.numeric(IDedTable[i, KMD_col]) - as.numeric(IDedTable[start, KMD_col]) > mz_win) {
      IDedTable[start:(i-1), KMD_col + 1] <- group
      group <- group + 1
      start <- i
    }
  }
  if (start != nrow_IDedTable) {
    IDedTable[start:nrow_IDedTable, KMD_col + 1] <- group
  }
  
  # Calculate homologous series
  series <- 1
  for (i in RowStartForFeatureTableData:(nrow_IDedTable - 1)) {
    if (is.na(IDedTable[i, KMD_col + 2]) || IDedTable[i, KMD_col + 2] != 0) {
      next
    }
    group <- as.numeric(IDedTable[i, KMD_col + 1])
    grow <- which(IDedTable[, KMD_col + 1] == group)
    mrow <- grow[which(((as.numeric(IDedTable[grow, KMD_col - 1]) - as.numeric(IDedTable[i, KMD_col - 1])) / Mass) %% 1 == 0)]
    members <- as.numeric(IDedTable[mrow, KMD_col - 1])
    if (length(unique(members)) < 2) {
      IDedTable[mrow, KMD_col + 2] <- NA
    } else {
      IDedTable[mrow, KMD_col + 2] <- series
      series <- series + 1
    }
  }
  IDedTable[which(is.na(IDedTable[, KMD_col + 2])), KMD_col + 2] <- 0
  
  # Sort by homologous series
  # IDedTable <- IDedTable[order(as.numeric(IDedTable[, ncol_IDedTable])),]
  
  m <- matrix(0, nrow(IDedTable), 2)
  IDedTable <- cbind(IDedTable, m)
  # Phantom columns
  ncol_IDedTable <- ncol(IDedTable) + 8
  colnames(IDedTable)[ncol_IDedTable - 9] <- paste(Names[x], "Homologous_series_count", sep="_")
  colnames(IDedTable)[ncol_IDedTable - 8] <- paste(Names[x], "MD_Filter", sep="_")
  hseries_count_cols <- append(hseries_count_cols, ncol_IDedTable - 9)
  
  # Mass defect filtering
  IDedTable[, ncol_IDedTable - 8] <- (as.numeric(IDedTable[, md_col]) >= lower & as.numeric(IDedTable[, md_col]) <= upper)
  
  # Pre-flagging sorting
  tmptable <- IDedTable[RowStartForFeatureTableData:nrow_IDedTable,]
  tmptable <- tmptable[order(as.numeric(tmptable[, ncol_IDedTable - 10]), as.numeric(tmptable[, ncol_IDedTable - 11]), as.numeric(tmptable[, ncol_IDedTable - 13]), as.numeric(tmptable[, Retention_col])),]
  IDedTable[RowStartForFeatureTableData:nrow_IDedTable,] <- tmptable
  
  # write.table(IDedTable, paste("C:/Users/pstel/Documents/PFAS Lab/Agilent/SMILES2MSMS/", "IDedTableIDed_FIN_KMD_scored_final_output_practice0.csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  
  # RT Flagging
  m <- matrix(0, nrow(IDedTable), 2)
  IDedTable <- cbind(IDedTable, m)
  ncol_IDedTable <- ncol_IDedTable - 1
  colnames(IDedTable)[ncol_IDedTable - 6] <- paste(Names[x], "RT_series", sep="_")
  series <- as.numeric(unique(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 9]))
  if (series[1] == 0) {
    series <- series[-1]
  }
  for (i in series) {
    curr <- which(IDedTable[, ncol_IDedTable - 9] == i)
    mseries <- as.numeric(unique(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]))
    if (length(mseries) < 2) {
      next
    }
    length_mseries <- length(mseries)
    for (j in 1:length_mseries) {
      mcurr <- curr[which(as.numeric(IDedTable[curr, ncol_IDedTable - 12]) == mseries[j])]
      if (mseries[j] == mseries[1]) {
        nex <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j + 1])]
        IDedTable[mcurr, ncol_IDedTable - 5] <- IDedTable[mcurr, Retention_col] < IDedTable[nex[length(nex)], Retention_col]
      } else if (mseries[j] == mseries[length(mseries)]) {
        prev <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j - 1])]
        IDedTable[mcurr, ncol_IDedTable - 6] <- IDedTable[mcurr, Retention_col] > IDedTable[prev[1], Retention_col]
      } else {
        prev <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j - 1])]
        IDedTable[mcurr, ncol_IDedTable - 6] <- IDedTable[mcurr, Retention_col] > IDedTable[prev[1], Retention_col]
        nex <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j + 1])]
        IDedTable[mcurr, ncol_IDedTable - 5] <- IDedTable[mcurr, Retention_col] < IDedTable[nex[length(nex)], Retention_col]
      }
    }
  }
  
  for (i in RowStartForFeatureTableData:nrow_IDedTable) {
    if (IDedTable[i, ncol_IDedTable - 6] == FALSE || IDedTable[i, ncol_IDedTable - 5] == FALSE) {
      IDedTable[i, ncol_IDedTable - 6] <- "RT not ordered"
    } else if (IDedTable[i, ncol_IDedTable - 6] == TRUE || IDedTable[i, ncol_IDedTable - 5] == TRUE) {
      IDedTable[i, ncol_IDedTable - 6] <- "RT ordered"
    } else {
      IDedTable[i, ncol_IDedTable - 6] <- "NA"
    }
  }
  
  IDedTable <- IDedTable[, 1:(ncol_IDedTable - 6)]
  
  ncol_IDedTable <- ncol_IDedTable + 1
  
  ## For debugging
  # write.table(IDedTable, file.path(OutDir, "NegIDed_test_optimized.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  
  if (RT_flagging == TRUE) {
    # Homologous seris count
    series <- as.numeric(unique(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 10]))
    if (series[1] == 0) {
      series <- series[-1]
    }
    none <- 2
    for (i in series) {
      curr <- which(IDedTable[, ncol_IDedTable - 10] == i)
      stable <- IDedTable[curr, c(ncol_IDedTable - 13, ncol_IDedTable - 7)]
      # Remove not ordered rows
      stable <- stable[!(is.na(stable[, 2]) == FALSE & stable[, 2] != "RT ordered" & stable[, 2] != "RT ordered, RT ordered"),]
      if (is.null(nrow(stable))) {
        none <- 1
      } else if (nrow(stable) == 0 && ncol(stable) == 2) {
        none <- 0
      }
      if (none == 2) {
        svec <- unique(stable[, 1])
        IDedTable[curr, ncol_IDedTable - 9] <- length(svec)
      } else if (none == 1) {
        IDedTable[curr, ncol_IDedTable - 9] <- 1
        none <- 2
      } else if (none == 0) {
        IDedTable[curr, ncol_IDedTable - 9] <- 0
        none <- 2
      }
    }
  } else {
    series <- as.numeric(unique(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 10]))
    if (series[1] == 0) {
      series <- series[-1]
    }
    for (i in series) {
      curr <- which(IDedTable[, ncol_IDedTable - 10] == i)
      svec <- unique(IDedTable[curr, ncol_IDedTable - 13])
      IDedTable[curr, ncol_IDedTable - 9] <- length(svec)
    }
  }
  
  # DT Flagging
  m <- matrix(0, nrow(IDedTable), 2)
  IDedTable <- cbind(IDedTable, m)
  ncol_IDedTable <- ncol_IDedTable - 1
  colnames(IDedTable)[ncol_IDedTable - 5] <- paste(Names[x], "DT_series", sep="_")
  series <- as.numeric(unique(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 9]))
  if (series[1] == 0) {
    series <- series[-1]
  }
  for (i in series) {
    curr <- which(IDedTable[, ncol_IDedTable - 9] == i)
    mseries <- as.numeric(unique(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]))
    if (length(mseries) < 2) {
      next
    }
    length_mseries <- length(mseries)
    for (j in 1:length_mseries) {
      mcurr <- curr[which(as.numeric(IDedTable[curr, ncol_IDedTable - 12]) == mseries[j])]
      if (mseries[j] == mseries[1]) {
        nex <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j + 1])]
        IDedTable[mcurr, ncol_IDedTable - 4] <- IDedTable[mcurr, DT_col] < IDedTable[nex[length(nex)], DT_col]
      } else if (mseries[j] == mseries[length(mseries)]) {
        prev <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j - 1])]
        IDedTable[mcurr, ncol_IDedTable - 5] <- IDedTable[mcurr, DT_col] > IDedTable[prev[1], DT_col]
      } else {
        prev <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j - 1])]
        IDedTable[mcurr, ncol_IDedTable - 5] <- IDedTable[mcurr, DT_col] > IDedTable[prev[1], DT_col]
        nex <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j + 1])]
        IDedTable[mcurr, ncol_IDedTable - 4] <- IDedTable[mcurr, DT_col] < IDedTable[nex[length(nex)], DT_col]
      }
    }
  }
  
  for (i in RowStartForFeatureTableData:nrow_IDedTable) {
    if (IDedTable[i, ncol_IDedTable - 5] == FALSE || IDedTable[i, ncol_IDedTable - 4] == FALSE) {
      IDedTable[i, ncol_IDedTable - 5] <- "DT not ordered"
    } else if (IDedTable[i, ncol_IDedTable - 5] == TRUE || IDedTable[i, ncol_IDedTable - 4] == TRUE) {
      IDedTable[i, ncol_IDedTable - 5] <- "DT ordered"
    } else {
      IDedTable[i, ncol_IDedTable - 5] <- "NA"
    }
  }
  
  IDedTable <- IDedTable[, 1:(ncol_IDedTable - 5)]
  
  ncol_IDedTable <- ncol_IDedTable + 2
  
  # write.table(IDedTable, paste("C:/Users/pstel/Documents/PFAS Lab/Agilent/SMILES2MSMS/", "IDedTableIDed_FIN_KMD_scored_final_output_practice1.csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
}

# timestamp()

# write.table(IDedTable, file.path(OutDir, "Test.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")

Potential_IDs_col <- NA
Frags_col <- NA

# if (Lipid == FALSE) {
m <- matrix(0, nrow(IDedTable), 8)
IDedTable <- cbind(IDedTable, m)
ncol_IDedTable <- ncol(IDedTable)
colnames(IDedTable)[ncol_IDedTable - 7] <- "Confident_ID"
colnames(IDedTable)[ncol_IDedTable - 6] <- "2+_homologous_series"
colnames(IDedTable)[ncol_IDedTable - 5] <- "3+_homologous_series"
colnames(IDedTable)[ncol_IDedTable - 4] <- "Tentative_ID"
colnames(IDedTable)[ncol_IDedTable - 3] <- "F_containing"
colnames(IDedTable)[ncol_IDedTable - 2] <- "Exact_mass_match"
colnames(IDedTable)[ncol_IDedTable - 1] <- "Score"
colnames(IDedTable)[ncol_IDedTable] <- "Score_Description"

# Scoring system
for (i in RowStartForFeatureTableData:nrow_IDedTable) {
  if (substr(IDedTable[i, ID_Ranked_col], 1, 2) == "1_" || substr(IDedTable[i, ID_Ranked_col], 1, 2) == "2_") {
    IDedTable[i, ncol_IDedTable - 7] <- TRUE
  } else {
    IDedTable[i, ncol_IDedTable - 7] <- FALSE
  }
  IDedTable[i, ncol_IDedTable - 6] <- FALSE
  for (hsc_col in hseries_count_cols){
    if (as.numeric(IDedTable[i, hsc_col]) >= 2) {
      IDedTable[i, ncol_IDedTable - 6] <- TRUE
    }
  }
  IDedTable[i, ncol_IDedTable - 5] <- FALSE
  for (hsc_col in hseries_count_cols){
    if (as.numeric(IDedTable[i, hsc_col]) >= 3) {
      IDedTable[i, ncol_IDedTable - 5] <- TRUE
    }
  }
  ##JPK Added this commenting as it was throwing an error and we don't need it... For potential IDs... MS/MS... Frag Counting...
  # if (is.na(IDedTable[i, Potential_IDs_col])) {
  #   IDedTable[i, ncol_IDedTable - 4] <- FALSE
  # } else {
  #   IDedTable[i, ncol_IDedTable - 4] <- TRUE
  # }
  # Fvec <- strsplit(IDedTable[i, Frags_col], "")
  # if (is.na(Fvec) == FALSE) {
  #   for (j in 1:length(Fvec[[1]])) {
  #     if (Fvec[[1]][j] == 'F') {
  #       IDedTable[i, ncol_IDedTable - 3] <- TRUE
  #       break
  #     }
  #   }
  # }
  if (IDedTable[i, ncol_IDedTable - 3] != TRUE) {
    IDedTable[i, ncol_IDedTable - 3] <- FALSE
  }
  if (substr(IDedTable[i, ID_Ranked_col], 1, 2) == "4_") {
    IDedTable[i, ncol_IDedTable - 2] <- TRUE
  } else {
    IDedTable[i, ncol_IDedTable - 2] <- FALSE
  }
}

# timestamp()

for (i in RowStartForFeatureTableData:nrow_IDedTable) {
  if (IDedTable[i, ncol_IDedTable - 7] == "TRUE" && IDedTable[i, ncol_IDedTable - 6] == "TRUE") {
    IDedTable[i, ncol_IDedTable - 1] <- "A"
    IDedTable[i, ncol_IDedTable] <- "Lvl 2 Schymanski: Confident ID (class specific dominant fragments observed along with exact mass) and 2+ in homologous series"
  } else if (IDedTable[i, ncol_IDedTable - 7] == "TRUE") {
    IDedTable[i, ncol_IDedTable - 1] <- "A-"
    IDedTable[i, ncol_IDedTable] <- "Lvl 2 Schymanski: Confident ID (class specific dominant fragments and exact mass)"
  } else if (IDedTable[i, ncol_IDedTable - 4] == "TRUE" && IDedTable[i, ncol_IDedTable - 6] == "TRUE") {
    IDedTable[i, ncol_IDedTable - 1] <- "B+"
    IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass) and 2+ in homologous series. Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
  } else if (IDedTable[i, ncol_IDedTable - 3] == "TRUE" && IDedTable[i, ncol_IDedTable - 6] == "TRUE") {
    IDedTable[i, ncol_IDedTable - 1] <- "B"
    IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS common fragment (F containing) and exact mass) and 2+ in homologous series. Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
  } else if (IDedTable[i, ncol_IDedTable - 4] == "TRUE" || IDedTable[i, ncol_IDedTable - 3] == "TRUE") {
    IDedTable[i, ncol_IDedTable - 1] <- "B-"
    IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass OR 1+ common PFAS fragment (F containing), and exact mass). Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
  } else if ((IDedTable[i, ncol_IDedTable - 9] == "TRUE" || IDedTable[i, ncol_IDedTable - 2] == "TRUE") && IDedTable[i, ncol_IDedTable - 5] == "TRUE") {
    IDedTable[i, ncol_IDedTable - 1] <- "D+"
    IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: No ID, possible PFAS:  Mass defect falling within -0.11 and 0.12 OR exact mass match, and 3+ within homologous series"
  } else if (IDedTable[i, ncol_IDedTable - 5] == "TRUE") {
    IDedTable[i, ncol_IDedTable - 1] <- "D"
    IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: No ID, possible PFAS: 3+ within homologous series"
  } else if (IDedTable[i, ncol_IDedTable - 9] == "TRUE" || IDedTable[i, ncol_IDedTable - 2] == "TRUE") {
    IDedTable[i, ncol_IDedTable - 1] <- "D-"
    IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: No ID, possible PFAS: Mass defect falling within -0.11 and 0.12 OR exact mass match"
  } else {
    IDedTable[i, ncol_IDedTable - 1] <- "E"
    IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: likely not PFAS"
  }
}

# timestamp()

# for (i in RowStartForFeatureTableData:nrow_IDedTable) {
if (ParallelComputing == TRUE) {
  
  IDedTable <- foreach (i = RowStartForFeatureTableData:nrow_IDedTable, .combine = rbind) %dopar% {
    #print(i)
    if (IDedTable[i, ncol_IDedTable - 1] == "D+" || IDedTable[i, ncol_IDedTable - 1] == "D" || IDedTable[i, ncol_IDedTable - 1] == "D-") {
      for (hsc_col in hseries_count_cols) {
        if (as.numeric(IDedTable[i, hsc_col]) > 1) { ## TEST
          tmptable <- IDedTable[order(as.numeric(IDedTable[, hsc_col - 1])),]
          hsc_rows <- which(tmptable[, hsc_col - 1] == IDedTable[i, hsc_col - 1])
          start <- hsc_rows[1]
          stop <- hsc_rows[length(hsc_rows)]
          if (length(grep("A", tmptable[start:stop, ncol_IDedTable - 1])) > 0) {
            IDedTable[i, ncol_IDedTable - 1] <- "C+"
            IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: 2+ in homologous series, and at least one confident PFAS identification within homologous series (A- or higher grade)"
          } else if (length(grep("B", tmptable[start:stop, ncol_IDedTable - 1])) > 0 && IDedTable[i, ncol_IDedTable - 1] != "C+" && IDedTable[i, ncol_IDedTable - 1] != "C") {
            if (IDedTable[i, hsc_col] >= 3) {
              IDedTable[i, ncol_IDedTable - 1] <- "C"
              IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: 3+ in homologous series, and at least one highly likely PFAS identified within homologous series (B- or higher grade)"
            } else if (IDedTable[i, hsc_col] >= 2) {
              IDedTable[i, ncol_IDedTable - 1] <- "C-"
              IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: homologous series 2+ within, and at least one highly likely PFAS identified within homologous series (B- or higher grade)"
            }
          }
        }
      }
    }
    return(IDedTable[i,])
  }
} else {
  for (i in RowStartForFeatureTableData:nrow_IDedTable) {
    #print(i)
    if (IDedTable[i, ncol_IDedTable - 1] == "D+" || IDedTable[i, ncol_IDedTable - 1] == "D" || IDedTable[i, ncol_IDedTable - 1] == "D-") {
      for (hsc_col in hseries_count_cols) {
        if (as.numeric(IDedTable[i, hsc_col]) > 1) { ## TEST
          tmptable <- IDedTable[order(as.numeric(IDedTable[, hsc_col - 1])),]
          hsc_rows <- which(tmptable[, hsc_col - 1] == IDedTable[i, hsc_col - 1])
          start <- hsc_rows[1]
          stop <- hsc_rows[length(hsc_rows)]
          if (length(grep("A", tmptable[start:stop, ncol_IDedTable - 1])) > 0) {
            IDedTable[i, ncol_IDedTable - 1] <- "C+"
            IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: 2+ in homologous series, and at least one confident PFAS identification within homologous series (A- or higher grade)"
          } else if (length(grep("B", tmptable[start:stop, ncol_IDedTable - 1])) > 0 && IDedTable[i, ncol_IDedTable - 1] != "C+" && IDedTable[i, ncol_IDedTable - 1] != "C") {
            if (IDedTable[i, hsc_col] >= 3) {
              IDedTable[i, ncol_IDedTable - 1] <- "C"
              IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: 3+ in homologous series, and at least one highly likely PFAS identified within homologous series (B- or higher grade)"
            } else if (IDedTable[i, hsc_col] >= 2) {
              IDedTable[i, ncol_IDedTable - 1] <- "C-"
              IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: homologous series 2+ within, and at least one highly likely PFAS identified within homologous series (B- or higher grade)"
            }
          }
        }
      }
    }
  }
}

# timestamp()

# Change Ds to D+s when in series with D+s
for (i in RowStartForFeatureTableData:nrow_IDedTable) {
  #print(i)
  if (IDedTable[i, ncol_IDedTable - 1] == "D") {
    D_rows <- which(IDedTable[, ncol_IDedTable - 11] == IDedTable[i, ncol_IDedTable - 11])
    start <- D_rows[1]
    stop <- D_rows[length(D_rows)]
    for (j in start:stop) {
      if (IDedTable[j, ncol_IDedTable - 1] == "D+") {
        IDedTable[i, ncol_IDedTable - 1] <- "D+"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: No ID, possible PFAS:  Mass defect falling within -0.11 and 0.12 OR exact mass match, and 3+ within homologous series"
        break
      }
    }
  }
}

# timestamp()

# write.table(IDedTable, paste("C:/Users/pstel/Documents/PFAS Lab/Agilent/SMILES2MSMS/", "IDedTableIDed_FIN_KMD_scored_final.csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")

# Next seven lines for debugging
# rm(list = ls())
# IDedTable_dir <- "C:/Users/pstel/Documents/PFAS Lab/Agilent/SMILES2MSMS/IDedTableIDed_FIN_KMD_scored_final.csv"
# IDedTable <- read.csv(IDedTable_dir,sep=",")
# IDedTable <- as.matrix(IDedTable)
# Retention_col <- 3
# nrow_IDedTable <- nrow(IDedTable)
# ncol_IDedTable <- ncol(IDedTable)

# Remove tentative and F containing columns
IDedTable <- IDedTable[, -((ncol_IDedTable - 4):(ncol_IDedTable - 3))]
ncol_IDedTable <- ncol(IDedTable)

######### Sorting

m <- matrix(0, nrow(IDedTable), 4)
IDedTable <- cbind(IDedTable, m)
ncol_IDedTable <- ncol(IDedTable) - 4

# Find max series for each row
for (i in RowStartForFeatureTableData:nrow_IDedTable) {
  IDedTable[i, ncol_IDedTable + 2] <- max(which(IDedTable[i, hseries_count_cols] == max(IDedTable[i, hseries_count_cols])))
}

# Max series series number and count
m <- matrix(0, nrow_IDedTable, 1)
IDedTable <- cbind(IDedTable, m)
for (i in RowStartForFeatureTableData:nrow_IDedTable) {
  IDedTable[i, ncol_IDedTable + 3] <- as.numeric(IDedTable[i, hseries_count_cols[as.numeric(IDedTable[i, ncol_IDedTable + 2])] - 1])
  IDedTable[i, ncol_IDedTable + 4] <- as.numeric(IDedTable[i, hseries_count_cols[as.numeric(IDedTable[i, ncol_IDedTable + 2])]])
}

# Sort by max series count, max series, and series number
IDedTable <- IDedTable[order(-as.numeric(IDedTable[, ncol_IDedTable + 4]), -as.numeric(IDedTable[, ncol_IDedTable + 2]), as.numeric(IDedTable[, ncol_IDedTable + 3])),]

# Create rank
rank <- 1
IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5] <- rank
for (i in (RowStartForFeatureTableData + 1):nrow(IDedTable)) {
  if (as.numeric(IDedTable[i, ncol_IDedTable + 4]) == 0) {
    rank <- rank + 1
  } else if (IDedTable[i, ncol_IDedTable + 3] != IDedTable[i - 1, ncol_IDedTable + 3]) {
    rank <- rank + 1
  }
  IDedTable[i, ncol_IDedTable + 5] <- rank
}

# Score
for (i in RowStartForFeatureTableData:nrow(IDedTable)) {
  if (IDedTable[i, ncol_IDedTable - 1] == "A" || IDedTable[i, ncol_IDedTable - 1] == "A-") {
    IDedTable[i, ncol_IDedTable + 1] <- 1
  } else if (IDedTable[i, ncol_IDedTable - 1] == "B+" || IDedTable[i, ncol_IDedTable - 1] == "B" || IDedTable[i, ncol_IDedTable - 1] == "B-") {
    IDedTable[i, ncol_IDedTable + 1] <- 2
  } else if (IDedTable[i, ncol_IDedTable - 1] == "E") {
    IDedTable[i, ncol_IDedTable + 1] <- 4
  } else {
    IDedTable[i, ncol_IDedTable + 1] <- 3
  }
}

# Sort by rank and score
IDedTable <- IDedTable[order(as.numeric(IDedTable[, ncol_IDedTable + 5]), as.numeric(IDedTable[, ncol_IDedTable + 1])),]

# Adjust rank
old <- IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5]
if (IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 1] == 1) {
  IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5] <- as.numeric(IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5]) - 10000000
} else if (IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 1] == 2) {
  IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5] <- as.numeric(IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5]) - 1000000
} else if (IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 1] == 4) {
  IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5] <- as.numeric(IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5]) + 1000000
}
for (i in (RowStartForFeatureTableData + 1):nrow(IDedTable)) {
  if (as.numeric(IDedTable[i, ncol_IDedTable + 1]) == 4) {
    old <- IDedTable[i, ncol_IDedTable + 5]
    IDedTable[i, ncol_IDedTable + 5] <- as.numeric(IDedTable[i, ncol_IDedTable + 5]) + 1000000
  } else if ((IDedTable[i, ncol_IDedTable + 5] != old)) {
    old <- IDedTable[i, ncol_IDedTable + 5]
    if (as.numeric(IDedTable[i, ncol_IDedTable + 1]) == 1) {
      IDedTable[i, ncol_IDedTable + 5] <- as.numeric(IDedTable[i, ncol_IDedTable + 5]) - 10000000
    } else if (as.numeric(IDedTable[i, ncol_IDedTable + 1]) == 2) {
      IDedTable[i, ncol_IDedTable + 5] <- as.numeric(IDedTable[i, ncol_IDedTable + 5]) - 1000000
    }
  } else {
    old <- IDedTable[i, ncol_IDedTable + 5]
    IDedTable[i, ncol_IDedTable + 5] <- as.numeric(IDedTable[i - 1, ncol_IDedTable + 5])
  }
}

# Sort by rank and m/z
IDedTable <- IDedTable[order(as.numeric(IDedTable[, ncol_IDedTable + 5]), as.numeric(IDedTable[, Mass_col])),]

## Is this right? I think so.
# Final aesthetics
ncol_IDedTable <- ncol_IDedTable - 1
for (i in RowStartForFeatureTableData:nrow_IDedTable) {
  if (IDedTable[i, ncol_IDedTable - 6] == "TRUE" || IDedTable[i, ncol_IDedTable - 6] == " TRUE") {
    IDedTable[i, ncol_IDedTable - 1] <- "TRUE"
  }
}
# }

#ensure all column names are consistent for PowerBI
colnames(IDedTable)[rowID_col]<-"row.ID"
colnames(IDedTable)[CCS_col]<-"CCS"
colnames(IDedTable)[DT_col]<-"DT"

# New formatting
m <- matrix("", nrow_IDedTable, 5)
n <- matrix("", nrow_IDedTable, 4)
IDedTable <- cbind(m, IDedTable[, Mass_col], IDedTable[, Retention_col], n, IDedTable[, -c(Mass_col, Retention_col)])
colnames(IDedTable)[1] <- "Score"
colnames(IDedTable)[2] <- "SeriesType_Identifier"
colnames(IDedTable)[3] <- "Name_or_Class"
colnames(IDedTable)[4] <- "Formula"
colnames(IDedTable)[5] <- "SMILES"
colnames(IDedTable)[6] <- "m/z"
colnames(IDedTable)[7] <- "Retention Time"
colnames(IDedTable)[8] <- "Adduct"
colnames(IDedTable)[9] <- "Unique"
colnames(IDedTable)[10] <- "Score_Description"
colnames(IDedTable)[11] <- "Needs_Validation"

ncol_IDedTable <- ncol(IDedTable)
ID_Ranked_col <- 9 + ID_Ranked_col
Potential_IDs_col <- 9 + Potential_IDs_col
Frags_col <- 9 + Frags_col

# if (Lipid == FALSE && TWeen_pos == FALSE) {
IDedTable[, 1] <- IDedTable[, ncol_IDedTable - 6]
IDedTable[, 2] <- paste(Names[as.numeric(IDedTable[, ncol_IDedTable - 3])], IDedTable[, ncol_IDedTable - 2], sep="_") ## Error??
IDedTable[, 10] <- IDedTable[, ncol_IDedTable - 5]

for (i in RowStartForFeatureTableData:nrow_IDedTable) {
  if (IDedTable[i, 1] == "A" || IDedTable[i, 1] == "A-") {
    primary <- strsplit(strsplit(IDedTable[i, ID_Ranked_col], " | ", TRUE)[[1]][1], ";", TRUE)[[1]]
    if (length(primary) > 1) {
      name <- strsplit(primary[1], "")[[1]]
      ## PS 5/14/2023
      # if (length(which(name == "-")) > 0) {
      #   IDedTable[i, 3] <- paste(name[3:(which(name == "-")[2] - 1)], collapse = "")
      # } else {
      # message(name)
      IDedTable[i, 3] <- paste(name[3:length(name)], collapse = "")
      # }
      ##
      IDedTable[i, 4] <- primary[2]
      IDedTable[i, 5] <- primary[3]
      IDedTable[i, 8] <- primary[4]
      ## PS 5/14/2023
      all <- strsplit(strsplit(IDedTable[i, ID_Ranked_col], " | ", TRUE)[[1]], ";", TRUE)
      u_entries <- c()
      length_all <- length(all)
      for (j in 1:length_all) {
        if (length(grep("\\[", all[[j]][length(all[[j]])])) > 0) {
          end <- strsplit(all[[j]][length(all[[j]])], "\\[")[[1]]
          all[[j]][length(all[[j]])] <- end[1]
          all[[j]] <- append(all[[j]], paste("[", end[2], sep=""))
        }
        u_entries <- append(u_entries, all[[j]][length(all[[j]]) - 1])
      }
      u_entries <- unique(u_entries)
      if (length(u_entries) == 1) {
        IDedTable[i, 9] <- "Yes"
      } else {
        IDedTable[i, 9] <- "No"
      }
      ##
    }
    # } else if (!is.na(IDedTable[i, Potential_IDs_col])) {
    #   primary <- strsplit(strsplit(IDedTable[i, Potential_IDs_col], "|", TRUE)[[1]][1], ";", TRUE)[[1]]
    #   if (length(primary) > 1) {
    #     name <- strsplit(primary[1], "")[[1]]
    #     ## PS 5/14/2023
    #     if (length(which(name == "-")) > 0) {
    #       IDedTable[i, 3] <- paste(name[1:(which(name == "-")[2] - 1)], collapse = "")
    #     } else {
    #       # message(name)
    #       IDedTable[i, 3] <- paste(name[1:length(name)], collapse = "")
    #     }
    #     ##
    #     IDedTable[i, 4] <- primary[2]
    #     IDedTable[i, 5] <- primary[3]
    #     IDedTable[i, 8] <- primary[4]
    #   }
    #   u_entries <- unique(strsplit(IDedTable[i, Potential_IDs_col], "|", TRUE)[[1]])
    #   if (length(u_entries) == 1) {
    #     IDedTable[i, 9] <- "Yes"
    #   } else {
    #     IDedTable[i, 9] <- "No"
    #   }
    # } else if (length(grep("DTXSID", IDedTable[i, Frags_col - 1])) > 0) {
    #   all <- strsplit(IDedTable[i, Frags_col - 1], "|", TRUE)[[1]]
    #   primary <- strsplit(all[grep("DTXSID", all)][1], ";", TRUE)[[1]]
    #   if (length(primary) > 1) {
    #     IDedTable[i, 3] <- primary[4]
    #     IDedTable[i, 4] <- primary[2]
    #     IDedTable[i, 5] <- primary[1]
    #     IDedTable[i, 8] <- "[M-H]-"
    #   }
    #   all <- strsplit(all, ";", TRUE)
    #   u_entries <- c()
    #   for (j in 1:length(all)) {
    #     u_entries <- append(u_entries, all[[j]][1])
    #   }
    #   u_entries <- unique(u_entries)
    #   if (length(u_entries) == 1) {
    #     IDedTable[i, 9] <- "Yes"
    #   } else {
    #     IDedTable[i, 9] <- "No"
    #   }
  } else if (substr(IDedTable[i, ID_Ranked_col], 1, 2) == "4_") {
    primary <- strsplit(strsplit(IDedTable[i, ID_Ranked_col], " | ", TRUE)[[1]][1], ";", TRUE)[[1]]
    # if (length(grep(";\\[", IDedTable[i, ID_Ranked_col])) == 0) {
    #   end <- strsplit(primary[length(primary)], "\\[")[[1]]
    #   primary[length(primary)] <- end[1]
    #   primary <- append(primary, paste("[", end[2], sep=""))
    # }
    ## PS 5/14/2023
    name <- strsplit(primary[1], "")[[1]]
    # if (length(which(name == "-")) > 0) {
    #   IDedTable[i, 3] <- paste(name[3:(which(name == "-")[2] - 1)], collapse = "")
    # } else {
    # message(name)
    IDedTable[i, 3] <- paste(name[3:length(name)], collapse = "")
    # }
    IDedTable[i, 4] <- primary[2]
    IDedTable[i, 5] <- primary[3]
    IDedTable[i, 8] <- primary[4]
    ##
    ## PS 5/14/2023
    all <- strsplit(strsplit(IDedTable[i, ID_Ranked_col], " | ", TRUE)[[1]], ";", TRUE)
    u_entries <- c()
    length_all <- length(all)
    for (j in 1:length_all) {
      if (length(grep("\\[", all[[j]][length(all[[j]])])) > 0) {
        end <- strsplit(all[[j]][length(all[[j]])], "\\[")[[1]]
        all[[j]][length(all[[j]])] <- end[1]
        all[[j]] <- append(all[[j]], paste("[", end[2], sep=""))
      }
      u_entries <- append(u_entries, all[[j]][length(all[[j]]) - 1])
    }
    u_entries <- unique(u_entries)
    if (length(u_entries) == 1) {
      IDedTable[i, 9] <- "Yes"
    } else {
      IDedTable[i, 9] <- "No"
    }
    ##
  }
  if (IDedTable[i, 1] != "A" & IDedTable[i, 1] != "A-") {
    IDedTable[i, 11] <- "Yes"
  } else {
    IDedTable[i, 11] <- "No"
  }
}

# } else if (TWeen_pos == TRUE) {
#   IDedTable[, 1] <- IDedTable[, ncol_IDedTable - 6]
#   IDedTable[, 2] <- paste(Names[as.numeric(IDedTable[, ncol_IDedTable - 3])], IDedTable[, ncol_IDedTable - 2], sep="_") ## Error??
#   IDedTable[, 10] <- IDedTable[, ncol_IDedTable - 5]
#   
#   for (i in RowStartForFeatureTableData:nrow_IDedTable) {
#     # Slight name modification
#     if (IDedTable[i, 1] == "A" || IDedTable[i, 1] == "A-") {
#       primary <- strsplit(strsplit(IDedTable[i, ID_Ranked_col], " | ", TRUE)[[1]][1], ";", TRUE)[[1]]
#       if (length(primary) > 1) {
#         name <- strsplit(primary[1], "")[[1]]
#         IDedTable[i, 3] <- paste(name[3:length(name)], collapse = "")
#         IDedTable[i, 4] <- primary[2]
#         IDedTable[i, 5] <- primary[3]
#         IDedTable[i, 8] <- primary[4]
#         IDedTable[i, 9] <- IDedTable[i, ID_Ranked_col + 4]
#       }
#       # Slight name modification
#     } else if (!is.na(IDedTable[i, Potential_IDs_col])) {
#       primary <- strsplit(strsplit(IDedTable[i, Potential_IDs_col], "|", TRUE)[[1]][1], ";", TRUE)[[1]]
#       if (length(primary) > 1) {
#         name <- strsplit(primary[1], "")[[1]]
#         IDedTable[i, 3] <- primary[1]
#         IDedTable[i, 4] <- primary[2]
#         IDedTable[i, 5] <- primary[3]
#         IDedTable[i, 8] <- primary[4]
#       }
#       u_entries <- unique(strsplit(IDedTable[i, Potential_IDs_col], "|", TRUE)[[1]])
#       if (length(u_entries) == 1) {
#         IDedTable[i, 9] <- "Yes"
#       } else {
#         IDedTable[i, 9] <- "No"
#       }
#       ### Name modification
#     } else if (length(grep("PEG", IDedTable[i, Frags_col - 1])) > 0) {
#       all <- strsplit(IDedTable[i, Frags_col - 1], "|", TRUE)[[1]]
#       # Take first PEG for name
#       # primary <- strsplit(all[grep("PEG", all)][1], ";", TRUE)[[1]]
#       primary <- all[grep("PEG", all)]
#       if (length(primary) > 1) {
#         IDedTable[i, 3] <- primary[1]
#         IDedTable[i, 4] <- "NA"
#         IDedTable[i, 5] <- "NA"
#         # [M+?]+ fill in ? from PEG name
#         IDedTable[i, 8] <- paste("[M+", strsplit(strsplit(primary[1], "_", TRUE)[[1]][2], "+", TRUE)[[1]][1], "]+", sep = "")
#       }
#       # Don't split, just check for unique
#       u_entries <- unique(primary)
#       if (length(u_entries) == 1) {
#         IDedTable[i, 9] <- "Yes"
#       } else {
#         IDedTable[i, 9] <- "No"
#       }
#       ### Name modification
#     } else if (substr(IDedTable[i, ID_Ranked_col], 1, 2) == "4_") {
#       primary <- strsplit(IDedTable[i, ID_Ranked_col], " | ", TRUE)[[1]]
#       name <- strsplit(primary[1], "")[[1]]
#       IDedTable[i, 3] <- paste(name[3:length(name)], collapse = "")
#       IDedTable[i, 4] <- "NA"
#       IDedTable[i, 5] <- "NA"
#       IDedTable[i, 8] <- "[M+?]+"
#       
#       u_entries <- unique(primary)
#       if (length(u_entries) == 1) {
#         IDedTable[i, 9] <- "Yes"
#       } else {
#         IDedTable[i, 9] <- "No"
#       }
#     }
#     if (IDedTable[i, 1] != "A" & IDedTable[i, 1] != "A-") {
#       IDedTable[i, 11] <- "Yes"
#     } else {
#       IDedTable[i, 11] <- "No"
#     }
#   }
# } else if (Lipid == TRUE) {
#   ################################## LIPIDMATCH SCORING SYSTEM
#   IDedTable[, 1] <- substr(IDedTable[, ID_Ranked_col], 1, 1)
#   IDedTable[, 2] <- paste(Names[1], IDedTable[, ncol_IDedTable - 3], sep="_")
#   IDedTable[, 3] <- IDedTable[, ID_Ranked_col + 2]
#   IDedTable[, 4] <- "CH2"
#   IDedTable[, 5] <- "CC"
#   
#   IDedTable[, 8] <- IDedTable[, ID_Ranked_col + 3]
#   IDedTable[which(is.na(IDedTable[, 8])), 8] <- 0
#   
#   IDedTable[, 9] <- IDedTable[, ID_Ranked_col + 4]
#   IDedTable[which(is.na(IDedTable[, 9])), 9] <- "No"
#   
#   IDedTable[which(IDedTable[, 1] == 1), 10] <- "Annotated using class based rules from standards including by fatty acyl chain constituents with ~5% false positive rate at the level of class and C:DB"
#   IDedTable[which(IDedTable[, 1] == 2), 10] <- "Annotated using class based rules from standards and DIA data with ~5% false positive rate at the level of class and C:DB"
#   IDedTable[which(IDedTable[, 1] == 3), 10] <- "Annotated using class based rules from standards and DDA data but only annotated by class (fatty acid composition not known)"
#   IDedTable[which(IDedTable[, 1] == 4), 10] <- "accurate mass match to database with very high false positive rate (>> 80%) so only use as a starting point for compound identification"
#   IDedTable[which(IDedTable[, 1] == 5), 10] <- "no match"
#   
#   IDedTable[, 11] <- "Yes"
# }

# timestamp()

IDedTable <- IDedTable[, 1:(ncol_IDedTable - 5)]

# Remove whitespace
IDedTable[IDedTable == ""] <- NA

# Remove "Na" from formulas in column 4
IDedTable[, 4] <- gsub("Na", "", IDedTable[, 4])

####################################################################################

# if (Lipid == FALSE) {
#   # Create a second output file with more info on B scores
#   IDedTable_B <- IDedTable
#   m <- matrix(0, nrow(IDedTable_B), 5)
#   IDedTable_B <- cbind(IDedTable_B, m)
#   nrow_IDedTable_B <- nrow(IDedTable_B)
#   ncol_IDedTable_B <- ncol(IDedTable_B)
#   colnames(IDedTable_B)[ncol_IDedTable_B - 4] <- "Check_Common_Fragment"
#   colnames(IDedTable_B)[ncol_IDedTable_B - 3] <- "Number_of_F_Fragments"
#   colnames(IDedTable_B)[ncol_IDedTable_B - 2] <- "Formula"
#   colnames(IDedTable_B)[ncol_IDedTable_B - 1] <- "MD"
#   colnames(IDedTable_B)[ncol_IDedTable_B] <- "Likely_PFAS"
#   
#   for (i in RowStartForFeatureTableData:nrow_IDedTable) {
#     count <- 0
#     frags_to_check <- c("C2F5", "C3F7", "SO2F", "SO3F", "CO2F", "C2F5O", "CF3O", "C2F3O2", "C3F7O", "C3O2F7", "CF3", "CH3FNSO2", "PO2F", "PO2F2", "SF5")
#     for (j in frags_to_check) {
#       if (length(grep(j, IDedTable_B[i, Frags_col])) > 0) {
#         count <- count + 1
#       }
#     }
#     IDedTable_B[i, ncol_IDedTable_B - 4] <- count
#     if (!is.na(IDedTable_B[i, Frags_col + 1])) {
#       IDedTable_B[i, ncol_IDedTable_B - 3] <- strsplit(IDedTable_B[i, Frags_col + 1], split="|", fixed=TRUE)[[1]][1]
#     }
#     if (is.na(IDedTable_B[i, 4])) {
#       IDedTable_B[i, ncol_IDedTable_B - 2] <- FALSE
#     } else {
#       IDedTable_B[i, ncol_IDedTable_B - 2] <- TRUE
#     }
#     if (as.numeric(IDedTable_B[i, Frags_col + 4]) > -0.25 && as.numeric(IDedTable_B[i, Frags_col + 4]) < 0.1) {
#       IDedTable_B[i, ncol_IDedTable_B - 1] <- TRUE
#     } else {
#       IDedTable_B[i, ncol_IDedTable_B - 1] <- FALSE
#     }
#     if (as.numeric(IDedTable_B[i, ncol_IDedTable_B - 4]) > 0 || as.numeric(IDedTable_B[i, ncol_IDedTable_B - 3]) > 3) {
#       IDedTable_B[i, ncol_IDedTable_B] <- TRUE
#     } else {
#       IDedTable_B[i, ncol_IDedTable_B] <- FALSE
#     }
#     
#     # Assign B--
#     if (IDedTable_B[i, ncol_IDedTable_B] == FALSE && (IDedTable_B[i, 1] == "B+" || IDedTable_B[i, 1] == "B" || IDedTable_B[i, 1] == "B-")) {
#       IDedTable_B[i, 1] <- "B--"
#       IDedTable_B[i, ncol_IDedTable_B - 6] <- "B--"
#       IDedTable_B[i, 10] <- "Lvl 5 Schymanski: Tentative ID, possible PFAS (1+ PFAS fragment (F containing) and exact mass). Note the fragment(s) observed may have other fragment formula assignments not containing F. Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
#       IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 5 Schymanski: Tentative ID, possible PFAS (1+ PFAS fragment (F containing) and exact mass). Note the fragment(s) observed may have other fragment formula assignments not containing F. Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
#     }
#     
#     # Correct B score descriptions that have SMILES
#     if (IDedTable_B[i, 1] == "B+" && !is.na(IDedTable_B[i, 5])) {
#       IDedTable_B[i, 10] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass) and 2+ in homologous series. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
#       IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass) and 2+ in homologous series. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
#     } else if (IDedTable_B[i, 1] == "B" && !is.na(IDedTable_B[i, 5])) {
#       IDedTable_B[i, 10] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS common fragment (F containing) and exact mass) and 2+ in homologous series. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
#       IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS common fragment (F containing) and exact mass) and 2+ in homologous series. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
#     } else if (IDedTable_B[i, 1] == "B-" && !is.na(IDedTable_B[i, 5])) {
#       IDedTable_B[i, 10] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass OR 1+ common PFAS fragment (F containing), and exact mass). Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
#       IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass OR 1+ common PFAS fragment (F containing), and exact mass). Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
#     } else if (IDedTable_B[i, 1] == "B--" && !is.na(IDedTable_B[i, 5])) {
#       IDedTable_B[i, 10] <- "Lvl 3 Schymanski: Tentative ID, possible PFAS (1+ PFAS fragment (F containing) and exact mass). Note the fragment(s) observed may have other fragment formula assignments not containing F. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
#       IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 3 Schymanski: Tentative ID, possible PFAS (1+ PFAS fragment (F containing) and exact mass). Note the fragment(s) observed may have other fragment formula assignments not containing F. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
#     }
#     
#   }
#   
#   ###################### Final Sorting ###########################
#   
#   m <- matrix(0, nrow(IDedTable_B), 4)
#   IDedTable_B <- cbind(IDedTable_B, m)
#   ncol_IDedTable_B <- ncol(IDedTable_B) - 4
#   hseries_count_cols <- hseries_count_cols + 9
#   Mass_col <- Mass_col + 4
#   
#   # Find max series for each row
#   for (i in RowStartForFeatureTableData:nrow_IDedTable_B) {
#     IDedTable_B[i, ncol_IDedTable_B + 2] <- max(which(IDedTable_B[i, hseries_count_cols] == max(IDedTable_B[i, hseries_count_cols])))
#   }
#   
#   # Max series series number and count
#   m <- matrix(0, nrow_IDedTable_B, 1)
#   IDedTable_B <- cbind(IDedTable_B, m)
#   for (i in RowStartForFeatureTableData:nrow_IDedTable_B) {
#     IDedTable_B[i, ncol_IDedTable_B + 3] <- as.numeric(IDedTable_B[i, hseries_count_cols[as.numeric(IDedTable_B[i, ncol_IDedTable_B + 2])] - 1])
#     IDedTable_B[i, ncol_IDedTable_B + 4] <- as.numeric(IDedTable_B[i, hseries_count_cols[as.numeric(IDedTable_B[i, ncol_IDedTable_B + 2])]])
#   }
#   
#   # Sort by max series count, max series, and series number
#   IDedTable_B <- IDedTable_B[order(-as.numeric(IDedTable_B[, ncol_IDedTable_B + 4]), -as.numeric(IDedTable_B[, ncol_IDedTable_B + 2]), as.numeric(IDedTable_B[, ncol_IDedTable_B + 3])),]
#   
#   # Create rank
#   rank <- 1
#   IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- rank
#   for (i in (RowStartForFeatureTableData + 1):nrow(IDedTable_B)) {
#     if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 4]) == 0) {
#       rank <- rank + 1
#     } else if (IDedTable_B[i, ncol_IDedTable_B + 3] != IDedTable_B[i - 1, ncol_IDedTable_B + 3]) {
#       rank <- rank + 1
#     }
#     IDedTable_B[i, ncol_IDedTable_B + 5] <- rank
#   }
#   
#   # Score
#   for (i in RowStartForFeatureTableData:nrow(IDedTable_B)) {
#     if (IDedTable_B[i, ncol_IDedTable_B - 6] == "A" || IDedTable_B[i, ncol_IDedTable_B - 6] == "A-") {
#       IDedTable_B[i, ncol_IDedTable_B + 1] <- 1
#     } else if (IDedTable_B[i, ncol_IDedTable_B - 6] == "B+" || IDedTable_B[i, ncol_IDedTable_B - 6] == "B" || IDedTable_B[i, ncol_IDedTable_B - 6] == "B-") {
#       IDedTable_B[i, ncol_IDedTable_B + 1] <- 2
#     } else if (IDedTable_B[i, ncol_IDedTable_B - 6] == "B--") {
#       IDedTable_B[i, ncol_IDedTable_B + 1] <- 3
#     } else if (IDedTable_B[i, ncol_IDedTable_B - 6] == "E") {
#       IDedTable_B[i, ncol_IDedTable_B + 1] <- 5
#     } else {
#       IDedTable_B[i, ncol_IDedTable_B + 1] <- 4
#     }
#   }
#   
#   # Sort by rank and score
#   IDedTable_B <- IDedTable_B[order(as.numeric(IDedTable_B[, ncol_IDedTable_B + 5]), as.numeric(IDedTable_B[, ncol_IDedTable_B + 1])),]
#   
#   # Adjust rank
#   old <- IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]
#   if (IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 1] == 1) {
#     IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]) - 10000000
#   } else if (IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 1] == 2) {
#     IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]) - 1000000
#   } else if (IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 1] == 3) {
#     IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]) - 100000
#   } else if (IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 1] == 5) {
#     IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]) + 1000000
#   }
#   for (i in (RowStartForFeatureTableData + 1):nrow(IDedTable_B)) {
#     if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 1]) == 5) {
#       old <- IDedTable_B[i, ncol_IDedTable_B + 5]
#       IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i, ncol_IDedTable_B + 5]) + 1000000
#     } else if ((IDedTable_B[i, ncol_IDedTable_B + 5] != old)) {
#       old <- IDedTable_B[i, ncol_IDedTable_B + 5]
#       if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 1]) == 1) {
#         IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i, ncol_IDedTable_B + 5]) - 10000000
#       } else if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 1]) == 2) {
#         IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i, ncol_IDedTable_B + 5]) - 1000000
#       } else if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 1]) == 3) {
#         IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i, ncol_IDedTable_B + 5]) - 100000
#       }
#     } else {
#       old <- IDedTable_B[i, ncol_IDedTable_B + 5]
#       IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i - 1, ncol_IDedTable_B + 5])
#     }
#   }
#   
#   # Sort by rank and m/z
#   IDedTable_B <- IDedTable_B[order(as.numeric(IDedTable_B[, ncol_IDedTable_B + 5]), as.numeric(IDedTable_B[, Mass_col])),]
#   
#   IDedTable_B <- IDedTable_B[, 1:ncol_IDedTable_B]
#   
#   IDedTable <- IDedTable_B
# } else if (Lipid == TRUE) {
#   m <- matrix(0, nrow(IDedTable), 1)
#   IDedTable <- cbind(IDedTable, m)
#   nrow_IDedTable <- nrow(IDedTable)
#   ncol_IDedTable <- ncol(IDedTable)
#   colnames(IDedTable)[ncol_IDedTable] <- "Number_of_F_Fragments"
#   Num_Frags_col <- which(colnames(IDedTable) == "Num_Frags")
#   if (length(Num_Frags_col) > 0) {
#     IDedTable[, ncol_IDedTable] <- substr(IDedTable[, Num_Frags_col], 1, 1)
#   }
# }

# # Add convenient columns
# m <- matrix("", nrow(IDedTable_B), 5)
# IDedTable_B <- cbind(IDedTable_B, m)
# ncol_IDedTable_B <- ncol(IDedTable_B)
# colnames(IDedTable_B)[(ncol_IDedTable_B - 4):ncol_IDedTable_B] <- c("Checked_Viz", "TRUE_Viz", "Class_Viz", "Cnumb_Viz", "Comment_Viz")
# 
# 
# if (NegPos == "Neg") {
#   # write.table(IDedTable, file.path(OutDir, "NegIDed_FIN_KMD_scored.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
#   write.table(IDedTable_B, file.path(OutDir, "NegIDed_FIN.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
# } else if (NegPos == "Pos") {
#   # write.table(IDedTable, file.path(OutDir, "PosIDed_FIN_KMD_scored.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
#   write.table(IDedTable_B, file.path(OutDir, "PosIDed_FIN.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
# } else if (NegPos == "Combined") {
#   write.table(IDedTable_B, file.path(OutDir, "CombinedIDed_FIN.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
# }

# Add convenient columns (can add any dummy column needed here as well as placeholders)
m <- matrix("", nrow(IDedTable), 6)
IDedTable <- cbind(IDedTable, m)
ncol_IDedTable <- ncol(IDedTable)
colnames(IDedTable)[(ncol_IDedTable - 5):ncol_IDedTable] <- c("Files","Checked_Viz", "TRUE_Viz", "Class_Viz", "Cnumb_Viz", "Comment_Viz")

##Pull out adducts
for (i in 1:nrow(IDedTable)) {
  Adducts_i<-unlist(gregexpr(pattern='|',IDedTable[i,3],fixed=TRUE))
  if (!is.na(Adducts_i[1])) {
    Adduct_Name<-substr(IDedTable[i,3],Adducts_i[1]+1,Adducts_i[2]-1)
    Adduct_Name<-gsub("\\(","[",Adduct_Name)
    Adduct_Name<-gsub("\\)","]",Adduct_Name)
    IDedTable[i,8]<-Adduct_Name
  }
}


FinalExport_IDedFin<-paste0(dirname(IDedTable_dir),"/",NegPos,"IDed_FIN.csv")
# OutDir_parts <- strsplit(IDedTable_dir, split = "_")[[1]]
# OutDir_parts <- paste(OutDir_parts[-length(OutDir_parts)], collapse = "_")
# OutDir_full <- paste0(OutDir_parts, "_FIN.csv")

write.table(IDedTable, FinalExport_IDedFin, sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA", qmethod="double")

EICcpp <- TRUE #Whether to use the C++ version of the EIC / MS1 generation code, seems to give the same results, faster but less well tested

if (EICcpp==TRUE) {
  #These are distributed by us (thanks Jonathan)
  #cannot installed online
  ## Michael is it necessary to detach these packages? They are required by RMassBank
  # if("Rcpp" %in% (.packages())){
  #  detach("package:mzR", unload=TRUE) 
  #  detach("package:RMassBank", unload=TRUE) 
  #  detach("package:Rdisop", unload=TRUE)
  #  detach("package:Rcpp", unload=TRUE) 
  # }
  library("EICim")
  library("MS1im")
}

if (!require("BiocManager", quietly = TRUE))install.packages("BiocManager", repos = "http://cran.us.r-project.org")

##Added formula prediction and other packages here

if("remotes" %in% rownames(installed.packages()) == FALSE) {install.packages("remotes", repos = "http://cran.us.r-project.org")}
library("remotes")

if("MassTools" %in% rownames(installed.packages()) == FALSE) {remotes::install_github("mjhelf/MassTools")}
library("MassTools")

#library(MetaboCoreUtils)
if("enviPat" %in% rownames(installed.packages()) == FALSE) {install.packages("enviPat", repos = "http://cran.us.r-project.org")}
library(enviPat, quietly=T)

if("RMassBank" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("RMassBank")}
library("RMassBank")

if("MetaboCoreUtils" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("MetaboCoreUtils")}
library("MetaboCoreUtils")

if("sqldf" %in% rownames(installed.packages()) == FALSE) {install.packages("sqldf", repos = "http://cran.us.r-project.org")}
library("sqldf")

if("RSQLite" %in% rownames(installed.packages()) == FALSE) {install.packages("RSQLite", repos = "http://cran.us.r-project.org")}
library("RSQLite")

if("agricolae" %in% rownames(installed.packages()) == FALSE) {install.packages("agricolae", repos = "http://cran.us.r-project.org")}
library("agricolae")

if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages("data.table", repos = "http://cran.us.r-project.org")}
library(data.table)
if("SearchTrees" %in% rownames(installed.packages()) == FALSE) {install.packages("SearchTrees", repos = "http://cran.us.r-project.org")}
library(SearchTrees)
if("comprehenr" %in% rownames(installed.packages()) == FALSE) {install.packages("comprehenr", repos = "http://cran.us.r-project.org")}
library(comprehenr)
if("mzR" %in% rownames(installed.packages()) == FALSE)  {
  BiocManager::install("mzR")}
library(mzR)

#######################EIC and MS1 R Scripts###################################
TargetEIC_Only <- TRUE
target_mzXML = dirname(IDedTable_dir)
OutputDirectory <- dirname(IDedTable_dir)

Target_mzXML<-list.files(target_mzXML, pattern = ".mzML", full.names = FALSE)
# if(TargetEIC_Only==TRUE) {
#   Target_mzXML<-grep(Target_mzXML, pattern = "Target|target|TARGET|Blank|blank|BLANK", invert=FALSE, value=TRUE)
# }

InputLibrary <- dirname(RepeatingUnits_dir)


if (EICcpp==TRUE) {
  source(paste(InputLibrary,"/Scripts/EIC_MS1_fns.R",sep=""))
  source(paste(InputLibrary,"/Scripts/genEIC.R",sep=""))
  source(paste(InputLibrary,"/Scripts/genIsoTable.R",sep=""))
  source(paste(InputLibrary,"/Scripts/MS1Spectragen.R",sep=""))
} else {
  print("not currently supported to have non-c++ version")
}

source(paste(InputLibrary,"/Scripts/IsotopePercentages.R",sep=""))
source(paste(InputLibrary,"/Scripts/Kaufmann.R",sep=""))
source(paste(InputLibrary,"/Scripts/Stats.R",sep=""))
source(paste(InputLibrary,"/Scripts/Manual_Review.R",sep=""))
source(paste(InputLibrary,"/Scripts/Formula_Prediction.R",sep=""))

#for FeatureID_Cols
#Make consistent especially DT and rowID
#mz, rt, DT, rowID, Formula
ISOstring<-"13C3;N;S;Cl2;18O;Br2"
#Tolerances, first is for selecting DT within RT tolerance and vice versa, second 
# High RT Tolerance goes with low DT tolerance to create an RT EIC.
# High DT Tolerance goes with low RT tolerance to create a DT Mobiligram.
GroupCSVDirectory<-""
if (EICcpp==TRUE) {
    arguments = construct_EM_arguments(
    PrecursorMassAccuracy = mz_win
    ,RT_Tolerances = c(0.1, 0.5)
    ,DT_Tolerances = c(0.1, 3)
    ,isIM = TRUE
    ,OutputDirectory = OutputDirectory
    # ,FeatureID_Cols = c(6,7,12,14,4) # mz, rt, id, dt, formula
    ,FeatureID_Cols = c(6,7,11+rowID_col,11+DT_col-1,4) # mz, rt, id, dt, formula
    ,GroupCSVDirectory = GroupCSVDirectory
    ,isostring = ISOstring
    ,isotable = paste(InputLibrary,"/Scripts/secondary_isotopes.csv",sep="")
    ,path_to_mzXML_Files = target_mzXML
    )
} else {
  arguments = construct_EM_arguments(
    PrecursorMassAccuracy = mz_win
    ,RT_Window = RT_Window
    ,OutputDirectory = OutputDirectory
    ,FeatureID_Cols = c(7,8,1,5)
    ,GroupCSVDirectory = GroupCSVDirectory
    ,isostring = ISOstring
    ,isotable = paste(InputLibrary,"/Scripts/secondary_isotopes.csv",sep="")
  )
  arguments$path_to_mzXML_Files = target_mzXML
}

if (NegPos=="Neg") {
  runallEM(Target_mzXML, arguments, runmode = TRUE, isNeg = TRUE)
} else if (NegPos=="Pos") {
  runallEM(Target_mzXML, arguments, runmode = TRUE, isNeg = FALSE)
} else {print("no recognized NegPos character string, should be Neg or Pos")}

#######################Kaufmann###################################
##Slightly edited from Modular, as RowStartForFeatureTableData is actually -1, so had to add 1, due to Paul's code
#Kaufmann eC, and ratios, values
MD_col_name <- "mass.defect"
MZ_col_name <- "m.z"
C13_col_name <- "13C1"

if (NegPos=="Pos") {
  outputFile<-DetermineFilePath("Pos")
  Kaufmann_eCs(outputFile,MD_col_name,MZ_col_name,C13_col_name)
  print("Kaufmann plot information complete, positive mode")
  #Column Headers for Editing and Markup
  Review_Col_Names(outputFile)
}
if (NegPos=="Neg") {
  outputFile<-DetermineFilePath("Neg")
  Kaufmann_eCs(outputFile,MD_col_name,MZ_col_name,C13_col_name)
  print("Kaufmann plot information complete, negative mode")
  #Column Headers for Editing and Markup
  Review_Col_Names(outputFile)
}

#######################Formula Prediction###################################
# GLOBALS
eList1 = c('C','H','N','O','S','F','Br','Cl')
eList2 = c('C','H','N','O','S','F','P')

# adducts <- data.frame(Adduct=c('[M-H]-','[M+H]+'),
#                       a=c(NA,'H1'),
#                       d=c('H1',NA))
adducts <- data.frame(Adduct=c('[M-H]-','[M+H]+','[M-H-CO2]-'),
                      a=c(NA,'H1',NA),
                      d=c('H1',NA,'C1O2H1'))

MF_topN=10  #Choose how many MFs to store
ppm_Window<-mz_ppm_win
# Formula Prediction can take a significant amount of time, if set to 0 the algorithm will try and predict formula for all potential PFAS (filtered by score)
Override_Predict<-1000 #set to 2 or higher to only predict formula for "n" top ranked features

if (NegPos=="Pos") {
  fh_Feature_MS1s<-paste0(OutputDirectory,"/Pos_Feature_MS1s.csv")
  fh_Feature_IDList<-paste0(OutputDirectory,"/PosIDed_FIN.csv")
  fh_MS1<-paste0(OutputDirectory,'/Pos_Feature_pFormula_MS1s.csv')
  fh_IDed_FIN<-paste0(OutputDirectory,'/PosIDed_FIN.csv')
  Formula_Prediction(Override_Predict,fh_Feature_MS1s,fh_Feature_IDList,fh_MS1,fh_IDed_FIN,MF_topN,adducts,eList1,eList2,q=1,Poltxt="+",(ppm_Window/2))
}
if (NegPos=="Neg") {
  fh_Feature_MS1s<-paste0(OutputDirectory,"/Neg_Feature_MS1s.csv")
  fh_Feature_IDList<-paste0(OutputDirectory,"/NegIDed_FIN.csv")
  fh_MS1<-paste0(OutputDirectory,'/Neg_Feature_pFormula_MS1s.csv')
  fh_IDed_FIN<-paste0(OutputDirectory,'/NegIDed_FIN.csv')
  Formula_Prediction(Override_Predict,fh_Feature_MS1s,fh_Feature_IDList,fh_MS1,fh_IDed_FIN,MF_topN,adducts,eList1,eList2,q=-1,Poltxt="-",ppmTol=(ppm_Window/2))
}
