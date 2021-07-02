library(ggbiplot)

# read in the dataframe (need to set row names)
df <- read.csv('chem_clean.csv', header = T, sep = ',', row.names = 1)
df

# remove Depth column
ns <- df[,-1]
head(ns)

# compute the principal components (centre and scale the data)
ns_pca <- prcomp(ns, center = T, scale. = T)
summary(ns_pca)
ns_pca$rotation
ns_pca$x

ns_site <- c(rep('Mt Hutt R. crithmifolius', 8), rep('Mt Hutt R. monroi', 8), rep("Porter's Pass", 8), rep('Lake Tennyson', 2))

ggbiplot(ns_pca, ellipse=TRUE, obs.scale = 1, var.scale = 1, groups=ns_site, varname.size = 3) +
  scale_colour_manual(name="Site", values= c("#b19f8e", "#83A552", "#5FA1F7", "#9B1F1A")) +
  theme_grey(base_size = 12) + # axes and legend font size 
  xlim(-2.8, 7.0)+ylim(-5.4, 4.5)




# remove porters_3 0-5 outlier and lobulatus samples
new <- ns[-c(21,25,26),]

new_pca <- prcomp(new, center = T, scale. = T)
summary(new_pca)
new_pca$rotation

new_site <- c(rep('Mt Hutt R. crithmifolius', 8), rep('Mt Hutt R. monroi', 8), rep("Porter's Pass", 7))

ggbiplot(new_pca, ellipse=TRUE, obs.scale = 1, var.scale = 1, groups=new_site, varname.size = 3, labels.size = 1) +
  scale_colour_manual(name="Site", values= c("#83A552", "#5FA1F7","#9B1F1A")) +
  theme_grey(base_size = 12) + 
  xlim(-3.2, 5.3)+ylim(-3.6, 3.6)







