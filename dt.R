
#Os dados são lidos como adquiridos por TAP, em formato csv
dpath = "~/astrodata/catalog/xdata.csv"
ncol = rep("numeric", 23)
csvDataset = read.csv(dpath, row.names = 'gaia_id', colClasses = c("character", ncol))
colours = c('FUV-r', 'NUV-r', 'u-r', 'J378-r', 'J395-r', 'J410-r', 'J430-r', 'g-r', 'J515-r', 'bp-r', 'r-J660', 'r-G','r-i','r-rp', 'r-J861', 'r-z', 'r-J', 'r-H', 'r-Ks', 'r-W1', 'r-W2', 'r-W3', 'r-W4')
colnames(csvDataset) = colours
label = row.names(csvDataset)
for (i in 1:length(label)) {label[i] = ifelse(substr(label[i],1,1) == "t", 2, 1)}  # 1 ssp, 2 star
label = as.integer(label)
dataset = cbind(csvDataset,label)

##Still need to rewrite the lines below

#dataset1 = data.frame(cbind(dataset$label,dataset$`u-r`,dataset$`F378-r`,dataset$`F395-r`,dataset$`F410-r`,dataset$`F430-r`,
                            dataset$`g-r`,dataset$`F515-r`,dataset$`r-F660`,dataset$`r-i`,dataset$`r-F861`,dataset$`r-z`))
#colnames(dataset1) = c('label','u-r','F378-r','F395-r','F410-r','F430-r','g-r','F515-r','r-F660','r-i','r-F861','r-z')
#row.names(dataset1) = row.names(dataset)

#write.table(dataset1, file="anndataset_n_l.txt")  # normalized, labeled
