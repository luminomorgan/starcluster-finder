# The data is acquired by TAP via Virtual Observatory software (TOPCAT) and will be reorganized so it can be properly used by the main algorithm.

dpath = "~/starcluster-finder/catalog/xdata.csv"
ncol = rep("numeric", 23)
csvRaw = read.csv(dpath, colClasses = c("character", ncol))
csvRaw = csvRaw[!duplicated(csvRaw$gaia_id), ]
dataset = csvRaw[,2:24]
row.names(dataset) = csvRaw$gaia_id
colours = c('FUV-r', 'NUV-r', 'u-r', 'J378-r', 'J395-r', 'J410-r', 'J430-r', 'g-r', 'J515-r', 'bp-r', 'r-J660', 'r-G','r-i','r-rp', 'r-J861', 'r-z', 'r-J', 'r-H', 'r-Ks', 'r-W1', 'r-W2', 'r-W3', 'r-W4')
colnames(dataset) = colours
label = row.names(dataset)
for (i in 1:length(label)) {label[i] = ifelse(substr(label[i],1,1) == "t", 2, 1)}  # 1 ssp, 2 star
label = as.integer(label)
dataset = cbind(dataset,label)

dataset1 = data.frame(cbind(dataset$label,dataset[1:23]))
colnames(dataset1) = c('label',colours)

write.table(dataset1, file = "~/starcluster-finder/catalog/xdata_nl.txt")  # normalized, labeled
