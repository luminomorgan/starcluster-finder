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
for (i in 1:length(label)) {label[i] = ifelse(substr(label[i],1,1) == "t", 1, 2)}  # 1 ssp, 2 star
label = as.integer(label)
dataset = cbind(dataset,label)

dataset1 = data.frame(cbind(dataset$label,dataset[1:23]))
colnames(dataset1) = c('label',colours)

sspPath = "~/starcluster-finder/synth/convolved_mags_04082020.dat"
sspData =  read.delim(sspPath, row.names = "X")
colnames(sspData) = c("H","J","Ks","J378","J395","J410","J430","J515","J660","J861","gSDSS","iSDSS","rSDSS","uSDSS","zSDSS","G","bp","rp","FUV","NUV","W1","W2","W3","W4")
sspData = sspData[!sspData$"J" == Inf, ] #drop inf values
sspCor = data.frame(sspData$FUV - sspData$rSDSS, sspData$NUV - sspData$rSDSS, sspData$uSDSS - sspData$rSDSS, sspData$J378 - sspData$rSDSS, sspData$J395 - sspData$rSDSS,  sspData$J410 - sspData$rSDSS,  sspData$J430 - sspData$rSDSS, sspData$gSDSS - sspData$rSDSS, sspData$J515 - sspData$rSDSS, sspData$bp - sspData$rSDSS, sspData$rSDSS - sspData$J660, sspData$rSDSS - sspData$G, sspData$rSDSS - sspData$iSDSS, sspData$rSDSS - sspData$rp, sspData$rSDSS - sspData$J861, sspData$rSDSS - sspData$zSDSS, sspData$rSDSS - sspData$J, sspData$rSDSS - sspData$H, sspData$rSDSS - sspData$Ks, sspData$rSDSS - sspData$W1, sspData$rSDSS - sspData$W2, sspData$rSDSS - sspData$W3, sspData$rSDSS - sspData$W4, row.names = row.names(sspData))
colnames(sspCor) = colours
ssplabel = row.names(sspCor)
for (i in 1:length(ssplabel)) {ssplabel[i] = ifelse(substr(ssplabel[i],1,1) == "t", 2, 1)}  # 1 ssp, 2 star
ssplabel = as.integer(ssplabel)
sspCor = cbind(ssplabel, sspCor)
colnames(sspCor) = c('label', colours)
df2 = rbind(dataset1, sspCor)

write.table(df2, file = "~/starcluster-finder/catalog/xdata+ssp.txt")  # normalized, labeled
