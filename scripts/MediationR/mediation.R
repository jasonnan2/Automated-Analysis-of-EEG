library(lavaan)
library(readxl)
library(Dict)
library(openxlsx)

wb <- createWorkbook()
all_data <- read_excel("A:\ClimateLD\finalTables\scalpTables.xlsx",sheet='scalpChoiceAlpha')
