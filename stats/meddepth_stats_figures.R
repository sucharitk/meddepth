######
## for the meditation rivalry experiment investigate the EEG band power correlates of subjective reports of meditation depth

library(car) # companion to applied regression
library(lme4) #Linear mixed effect model
library(emmeans) # Least squares means
library(ggpubr)
library(tidyverse)
library(gridExtra)
library(lmerTest)
library(ggplot2)
library(sjPlot) # for plot_model

### read the data file and initialise the factors
medriv_data = read.csv(
  '~/OneDrive/Projects/Experiments/Meditation_Rivalry/Data/icacomprem_med_medriv_physio_freqs_fft.csv',
  header = TRUE, sep = ',')
chanlabs <- read.csv(
  '~/OneDrive/Projects/Experiments/Meditation_Rivalry/Data/channames.csv',
  sep = ',')
channames = chanlabs$chans

# assign factors as factors
medriv_data$subj <- as.factor(medriv_data$subj)
medriv_data$group <- as.factor(medriv_data$group)

medriv_data$group <- as.character.factor(medriv_data$group)
medriv_data$group[medriv_data$group==1] <- "Meditators"
medriv_data$group[medriv_data$group==2] <- "Controls"
medriv_data$block <- as.factor(medriv_data$block)
medriv_data$subjcode <- as.factor(medriv_data$subjcode)
medriv_data$block <- as.character.factor(medriv_data$block)
medriv_data$block[medriv_data$block==1] <- "BL"
medriv_data$block[medriv_data$block==2] <- "CH"
medriv_data$block[medriv_data$block==3] <- "M1"
medriv_data$block[medriv_data$block==4] <- "M2"

nchans = max(medriv_data$chans)

medriv_data$hrate[medriv_data$subjcode==4160] = NA
medriv_data$hrv[medriv_data$subjcode==4160] = NA
medriv_data$hrate[medriv_data$subjcode==4158 & (medriv_data$block=="M2") ] = NA
medriv_data$hrv[medriv_data$subjcode==4158 & (medriv_data$block=="M2") ] = NA

medriv_data$rrate[medriv_data$iqrr>5] = NA

stat_summ_size <- .3
font_size <- 10
### plot of self-reports of 5 dimensions of MEDEQ
#####
qclusts = c(1:5)
nqclusts = length(qclusts)
gp.ph <- list()
tstnum = matrix(nrow=2, ncol=nqclusts)
pltlab = c("A", "B", "C", "D", "E")
ylabs = c("DL0: Hindrances", "DL1: Relaxation", 
          "DL2: Concentration", "DL3: Essential qualities", "DL4: Nonduality")

for (qc in 1:nqclusts)
{
  cname = paste("q",qc,sep = "")
  md.plot <- medriv_data %>% 
    dplyr::select(cname,"group","subj","block") %>%
    group_by(group, subj, block) %>%
    dplyr::rename(qr = cname) %>%
    dplyr::summarise(qr = mean(qr, na.rm = TRUE))

    gp.ph[[qc]] <- ggplot(md.plot, aes(x = block, y = qr
                                     ,colour = block # comment this line for greyscale
                                     )) +
    facet_wrap(~group) +
    theme_pubclean() +
    stat_summary(fun.data = mean_se, size = stat_summ_size) +
    coord_cartesian(ylim = c(1.8, 11.5)) + 
    scale_y_continuous(breaks=seq(2,12,by=2)) +
    xlab("Block") + ylab(ylabs[qc]) +
    theme(legend.position = "none", 
          text = element_text(size=font_size), 
          plot.caption = element_text(size = font_size),
          plot.tag = element_text(size = font_size, face = "bold")) +
    labs(tag = pltlab[qc]) 
    # + scale_color_grey() # uncomment this line for greyscale
}

astsize <- 3.5
nssize <- 1.8
ycomp <- c(8.3, 9.5, 10.7)
ysigast <- c(8.4, 9.6, 10.8)
ysigns <- c(8.6, 9.8, 11)

labtext1 = data.frame(  label = c("*", "*"),  group = c("Controls", "Meditators"))
my_comparisons <- list(c("BL", "M1"))
gp.ph[[2]] <- gp.ph[[2]] + stat_compare_means(comparisons = my_comparisons, 
                              label.y = ycomp[2], method = "wilcox.test",
                              na.rm = TRUE, paired = TRUE,
                              symnum.args = symnum.args, tip.length = .0) +
  geom_text(data = labtext1, mapping =aes(x=2, y=ysigast[2], label=label),
            size = astsize, colour = "black") 

labtext1 = data.frame(  label = c("*", "*"),  group = c("Controls", "Meditators"))
labtext2 = data.frame(  label = c("****", "****"),  group = c("Controls", "Meditators"))
labtext3 = data.frame(  label = c("****", "****"),  group = c("Controls", "Meditators"))
my_comparisons <- list( c("BL", "CH"), c("BL", "M1"), c("BL", "M2"))
gp.ph[[3]] <- gp.ph[[3]] + 
  stat_compare_means(comparisons = my_comparisons,
                     label.y = ycomp, method = "wilcox.test",
                     na.rm = TRUE, paired = TRUE,
                     symnum.args = symnum.args, tip.length = .0) +
  geom_text(data = labtext1, mapping =aes(x=1.5, y=ysigast[1], label=label),
            size = astsize, colour = "black") +
  geom_text(data = labtext2, mapping =aes(x=2, y=ysigast[2], label=label),
            size = astsize, colour = "black") +
  geom_text(data = labtext2, mapping =aes(x=2.5, y=ysigast[3], label=label),
            size = astsize, colour = "black") 

labtext1 = data.frame(  label = c("", "****"),  group = c("Controls", "Meditators"))
labtext2 = data.frame(  label = c("ns", ""),  group = c("Controls", "Meditators"))
labtext3 = data.frame(  label = c("", "****"),  group = c("Controls", "Meditators"))
labtext4 = data.frame(  label = c("ns", ""),  group = c("Controls", "Meditators"))
labtext5 = data.frame(  label = c("", "****"),  group = c("Controls", "Meditators"))
labtext6 = data.frame(  label = c("ns", ""),  group = c("Controls", "Meditators"))
gp.ph[[4]] <- gp.ph[[4]] + 
  stat_compare_means(comparisons = my_comparisons,
                     label.y = ycomp, method = "wilcox.test",
                     na.rm = TRUE, paired = TRUE,
                     symnum.args = symnum.args, tip.length = .0) +
  geom_text(data = labtext1, mapping =aes(x=1.5, y=ysigast[1], label=label),
            size = astsize, colour = "black") +
  geom_text(data = labtext2, mapping =aes(x=1.5, y=ysigns[1], label=label),
            size = nssize, colour = "black") + 
  geom_text(data = labtext3, mapping =aes(x=2, y=ysigast[2], label=label),
            size = astsize, colour = "black") +
  geom_text(data = labtext4, mapping =aes(x=2, y=ysigns[2], label=label),
            size = nssize, colour = "black") + 
  geom_text(data = labtext3, mapping =aes(x=2.5, y=ysigast[3], label=label),
            size = astsize, colour = "black") +
  geom_text(data = labtext4, mapping =aes(x=2.5, y=ysigns[3], label=label),
            size = nssize, colour = "black")

labtext1 = data.frame(  label = c("", "****"),  group = c("Controls", "Meditators"))
labtext2 = data.frame(  label = c("ns", ""),  group = c("Controls", "Meditators"))
labtext3 = data.frame(  label = c("", "****"),  group = c("Controls", "Meditators"))
labtext4 = data.frame(  label = c("ns", ""),  group = c("Controls", "Meditators"))
labtext5 = data.frame(  label = c("", "****"),  group = c("Controls", "Meditators"))
labtext6 = data.frame(  label = c("ns", ""),  group = c("Controls", "Meditators"))
gp.ph[[5]] <- gp.ph[[5]] + 
  stat_compare_means(comparisons = my_comparisons,
                     label.y = ycomp, method = "wilcox.test",
                     na.rm = TRUE, paired = TRUE,
                     symnum.args = symnum.args, tip.length = .0) +
  geom_text(data = labtext1, mapping =aes(x=1.5, y=ysigast[1], label=label),
            size = astsize, colour = "black") +
  geom_text(data = labtext2, mapping =aes(x=1.5, y=ysigns[1], label=label),
            size = nssize, colour = "black") + 
  geom_text(data = labtext3, mapping =aes(x=2, y=ysigast[2], label=label),
            size = astsize, colour = "black") +
  geom_text(data = labtext4, mapping =aes(x=2, y=ysigns[2], label=label),
            size = nssize, colour = "black") + 
  geom_text(data = labtext3, mapping =aes(x=2.5, y=ysigast[3], label=label),
            size = astsize, colour = "black") +
  geom_text(data = labtext4, mapping =aes(x=2.5, y=ysigns[3], label=label),
            size = nssize, colour = "black")

### figure 1
grid.arrange(grobs = c(gp.ph), ncol = 3) # fig size 6.85 x 4.5 inch



# statistical test for the 5 meditation depths
md.plot <- medriv_data %>%
  group_by(group, subj, block) %>%
  summarise(q1 = mean(q1, na.rm = TRUE),
            q2 = mean(q2, na.rm = TRUE),
            q3 = mean(q3, na.rm = TRUE),
            q4 = mean(q4, na.rm = TRUE),
            q5 = mean(q5, na.rm = TRUE)) %>%
  pivot_longer(
    cols = starts_with("q"),
    names_to = "qnum",
    names_prefix = "q",
    values_to = "qresp"
  )
md.aov <- lmer(qresp ~ group*block*qnum + (1|subj), data = md.plot)
Anova(md.aov)

# statistical test for the 5 meditation depths
md.plot <- medriv_data %>%
  group_by(group, subj, block) %>%
  summarise(qr = mean(q1, na.rm = TRUE))
md.aov <- lmer(qr ~ group*block + (1|subj), data = md.plot)
plot(md.aov)
qqPlot(resid(md.aov))
Anova(md.aov)
md.aov <- lmer(qr ~ group+block + (1|subj), data = md.plot)
plot(md.aov, abline = TRUE)
qqPlot(resid(md.aov))
Anova(md.aov)


# relaxation
md.plot <- filter(medriv_data) %>%
  group_by(group, subj, block) %>%
  summarise(qr = mean(q2, na.rm = TRUE))

md.aov <- lmer(qr ~ group*block + (1|subj), data = md.plot)
plot(md.aov)
qqPlot(resid(md.aov))
Anova(md.aov)
md.aov <- lmer(qr ~ group+block + (1|subj), data = md.plot)
plot(md.aov)
qqPlot(resid(md.aov))
Anova(md.aov)
lsmeans(md.aov, pairwise ~ block)$contrasts

# concentration
md.plot <- filter(medriv_data) %>%
  group_by(group, subj, block) %>%
  summarise(qr = mean(q3, na.rm = TRUE))
md.aov <- lmer(qr ~ group*block + (1|subj), data = md.plot)
plot(md.aov)
qqPlot(resid(md.aov))
Anova(md.aov)
md.aov <- lmer(qr ~ group+block + (1|subj), data = md.plot)
plot(md.aov)
qqPlot(resid(md.aov))
Anova(md.aov)
lsmeans(md.aov, pairwise ~ block)$contrasts
lsmeans(md.aov, pairwise ~ group)$contrasts

# transpersonal qualities
md.plot <- filter(medriv_data) %>%
  group_by(group, subj, block) %>%
  summarise(qr = mean(q4, na.rm = TRUE))
# md.plot <- na.omit(md.plot)
md.aov <- lmer(qr ~ group*block + (1|subj), data = md.plot)
plot(md.aov)
qqPlot(resid(md.aov))
Anova(md.aov)
# anova(md.aov)
lsmeans(md.aov, pairwise ~ block|group)$contrasts
lsmeans(md.aov, pairwise ~ group|block)$contrasts

# nonduality
md.plot <- filter(medriv_data) %>%
  group_by(group, subj, block) %>%
  summarise(qr = mean(q5, na.rm = TRUE))
md.aov <- lmer(qr ~ group*block + (1|subj), data = md.plot)
plot(md.aov)
qqPlot(resid(md.aov))
Anova(md.aov)
lsmeans(md.aov, pairwise ~ block|group)$contrasts
lsmeans(md.aov, pairwise ~ group|block)$contrasts



#####  
### test for full model block x group x qnum

lasttime = 15*60

md.plot <- filter(medriv_data
                  ,chans %in% which(!(channames %in% c("FT9", "FT10", "TP9", "TP10", "FP1", "Fp2")))
                  , block=="BL" | block=="CH" | timefrend2 < lasttime 
                  # block is either baseline or chanting (~7 min blocks) or 
                  # take the last 7 minutes of the meditation blocks
) %>% 
  group_by(group, subj, block, subjcode) %>%
  summarise(alpha = mean(alpha, na.rm = TRUE),
            beta = mean(beta, na.rm = TRUE),
            theta = mean(theta, na.rm = TRUE),
            gamma = mean(gamma, na.rm = TRUE),
            hrv = mean(hrv, na.rm = TRUE)/10,
            hrate = mean(hrate, na.rm = TRUE)/10,
            rrate = mean(rrate, na.rm = TRUE)/10,
            q1 = mean(q1, na.rm = TRUE),
            q2 = mean(q2, na.rm = TRUE),
            q3 = mean(q3, na.rm = TRUE),
            q4 = mean(q4, na.rm = TRUE),
            q5 = mean(q5, na.rm = TRUE))

# put all the questionnaires in a single column
md.plot <- md.plot %>%
  pivot_longer(
    cols = starts_with("q"),
    names_to = "qnum",
    names_prefix = "q",
    values_to = "qresp"
  )

# first evaluate the full model of the 4 freq-bands and 3 physio measures with the 3-way interactions 

chsq <- matrix(nrow = 13, ncol = 14)
inum <- 1:14


# first evaluate the full model of the 4 freq-bands and 3 physio measures with the 3-way interactions 
md.aov <- lmer(qresp ~ theta*qnum*group+alpha*qnum*group+beta*qnum*group+
                 gamma*qnum*group+
                 hrate*qnum*group+hrv*qnum*group+rrate*qnum*group
               + (1|subj) , md.plot, REML = FALSE)
# first evaluate the full model of the 4 freq-bands and 3 physio measures with the 3-way interactions 

aov <- Anova(md.aov)
chsq[1,inum] <- c(aov$Chisq[25:31], aov$Chisq[10], aov$Chisq[13], aov$Chisq[15]
              , aov$Chisq[17], aov$Chisq[19], aov$Chisq[21], aov$Chisq[23])
chsq[1,inum] <- c(aov$Chisq[25:31], aov$Chisq[10], aov$Chisq[13], aov$Chisq[15]
                  , aov$Chisq[17], aov$Chisq[19], aov$Chisq[21], aov$Chisq[23])
plot(md.aov)

# start by removing the 3-way ints and the associated 2-way ints with group
m2 <- update(md.aov, ~.-qnum:group:rrate)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-rrate:group)
md.aov <- m2 # accept the reduced model
Anova(md.aov)
aov <- Anova(md.aov)
inum <- inum[-7]
chsq[2,inum] <- c(aov$Chisq[24:29], aov$Chisq[10], aov$Chisq[13], aov$Chisq[15]
              , aov$Chisq[17], aov$Chisq[19], aov$Chisq[21], aov$Chisq[23])



m2 <- update(md.aov, ~.-qnum:group:hrate)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-hrate:group)
md.aov <- m2 # accept the reduced model
Anova(md.aov)
aov <- Anova(md.aov)
inum <- inum[-5]
chsq[3,inum] <- c(aov$Chisq[23:27], aov$Chisq[10], aov$Chisq[13], aov$Chisq[15]
                  , aov$Chisq[17], aov$Chisq[19], aov$Chisq[20], aov$Chisq[22])



m2 <- update(md.aov, ~.-qnum:group:theta)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-theta:group)
md.aov <- m2 # accept the reduced model
Anova(md.aov)
aov <- Anova(md.aov)
inum <- inum[-1]
chsq[4,inum] <- c(aov$Chisq[22:25], aov$Chisq[10], aov$Chisq[12], aov$Chisq[14]
                  , aov$Chisq[16], aov$Chisq[18], aov$Chisq[19], aov$Chisq[21])


m2 <- update(md.aov, ~.-qnum:group:alpha)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-alpha:group)
md.aov <- m2 # accept the reduced model
Anova(md.aov)
aov <- Anova(md.aov)
inum <- inum[-1]
chsq[5,inum] <- c(aov$Chisq[21:23], aov$Chisq[10], aov$Chisq[12], aov$Chisq[13]
                  , aov$Chisq[15], aov$Chisq[17], aov$Chisq[18], aov$Chisq[20])


m2 <- update(md.aov, ~.-qnum:group:gamma)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-gamma:group)
md.aov <- m2 # accept the reduced model
Anova(md.aov)
aov <- Anova(md.aov)
inum <- inum[-2]
chsq[6,inum] <- c(aov$Chisq[20:21], aov$Chisq[10], aov$Chisq[12], aov$Chisq[13]
                  , aov$Chisq[15], aov$Chisq[16], aov$Chisq[17], aov$Chisq[19])


m2 <- update(md.aov, ~.-qnum:group:beta)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-beta:group)
md.aov <- m2 # accept the reduced model
Anova(md.aov)
aov <- Anova(md.aov)
inum <- inum[-1]
chsq[7,inum] <- c(aov$Chisq[19], aov$Chisq[10], aov$Chisq[12], aov$Chisq[13]
                  , aov$Chisq[14], aov$Chisq[15], aov$Chisq[16], aov$Chisq[18])

m2 <- update(md.aov, ~.-qnum:group:hrv)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-hrv:group)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-qnum:group) # remove the last remaining group interaction with qnum
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-group) # remove group main
md.aov <- m2 # accept the reduced model
Anova(md.aov)
aov <- Anova(md.aov)
inum <- inum[-1]
chsq[8,inum] <- c(aov$Chisq[9:15])


## now remove 2-way ints
m2 <- update(md.aov, ~.-qnum:hrv)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-hrv)
md.aov <- m2 # accept the reduced model
Anova(md.aov)
aov <- Anova(md.aov)
inum <- inum[-6]
chsq[9,inum] <- c(aov$Chisq[8:13])

m2 <- update(md.aov, ~.-qnum:beta)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-beta)
md.aov <- m2 # accept the reduced model
Anova(md.aov)
aov <- Anova(md.aov)
inum <- inum[-3]
chsq[10,inum] <- c(aov$Chisq[7:11])

m2 <- update(md.aov, ~.-qnum:hrate)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-hrate)
md.aov <- m2 # accept the reduced model
Anova(md.aov)
aov <- Anova(md.aov)
inum <- inum[-4]
chsq[11,inum] <- c(aov$Chisq[6:9])

m2 <- update(md.aov, ~.-qnum:gamma)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-gamma)
md.aov <- m2 # accept the reduced model
Anova(md.aov)
aov <- Anova(md.aov)
inum <- inum[-3]
chsq[12,inum] <- c(aov$Chisq[5:7])

m2 <- update(md.aov, ~.-qnum:rrate)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-rrate)
md.aov <- m2 # accept the reduced model
Anova(md.aov)
aov <- Anova(md.aov)
inum <- inum[-3]
chsq[13,inum] <- c(aov$Chisq[4:5])

m2 <- update(md.aov, ~.-qnum:theta)
anova(md.aov,m2)
# do not accept the reduced model

m2 <- update(md.aov, ~.-qnum:alpha)
anova(md.aov,m2)
# do not accept the reduced model


#### final model

### this is the final model <- evaluate with reml = true to get accurate p-values
md.aov <- lmer(qresp ~ alpha*qnum + 
                 theta*qnum +
                 (1|subj) , md.plot, REML = TRUE)

plot(md.aov)
qqPlot(resid(md.aov))
Anova(md.aov)
emt = emtrends(md.aov, "qnum", var = "alpha")
emt    # estimated slopes
pairs(emt)    # pairwise comparisons
emt = emtrends(md.aov, "qnum", var = "theta")
emt    # estimated slopes
pairs(emt)    # pairwise comparisons


# Heatmap 
tmod <- (sqrt(chsq)) %>%
  as_tibble() %>%
  rowid_to_column(var="Iteration") %>%
  gather(key="Interaction", value="Chi", -1) %>%
  mutate(Interaction=as.numeric(gsub("V","",Interaction)))

tmod$Chi[is.na(tmod$Chi)] <- 0

tmod$Interaction = as.character(tmod$Interaction)
tmod$Interaction[tmod$Interaction=="1"] = "3-theta"
tmod$Interaction[tmod$Interaction=="2"] = "3-alpha"
tmod$Interaction[tmod$Interaction=="3"] = "3-beta"
tmod$Interaction[tmod$Interaction=="4"] = "3-gamma"
tmod$Interaction[tmod$Interaction=="5"] = "3-hr"
tmod$Interaction[tmod$Interaction=="6"] = "3-hrv"
tmod$Interaction[tmod$Interaction=="7"] = "3-rr"
tmod$Interaction[tmod$Interaction=="8"] = "2-theta"
tmod$Interaction[tmod$Interaction=="9"] = "2-alpha"
tmod$Interaction[tmod$Interaction=="10"] = "2-beta"
tmod$Interaction[tmod$Interaction=="11"] = "2-gamma"
tmod$Interaction[tmod$Interaction=="12"] = "2-hr"
tmod$Interaction[tmod$Interaction=="13"] = "2-hrv"
tmod$Interaction[tmod$Interaction=="14"] = "2-rr"

modred <- ggplot(tmod, aes(Iteration, Interaction, fill= Chi)) + 
  geom_tile() + 
  theme_pubclean() + 
  scale_y_discrete(limits=c("3-rr",
                            "3-hr",
                            "3-theta",
                            "3-alpha",
                            "3-gamma",
                            "3-beta",
                            "3-hrv",
                            "2-hrv",
                            "2-beta",
                            "2-hr",
                            "2-gamma",
                            "2-rr",
                            "2-theta",
                            "2-alpha")) + 
  scale_x_discrete(limits = c(1:13)) + labs(tag = "A") +
  theme(axis.text.y = element_text(hjust = 0),
        text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = font_size, face = "bold")) +
  xlab("Model Reduction Step #") +
  ylab("Interaction Type")
        


md2 <- md.plot %>%
  filter(group=="Meditators") %>%
  dplyr::select(qresp, alpha, theta, qnum)
md.aov <- lmer(qresp ~ alpha*qnum +
                 theta*qnum + 
                 (1|subj) , md2)
plot(md.aov)
qqPlot(resid(md.aov))
Anova(md.aov)

md2 <- md.plot %>%
  filter(group=="Controls") %>%
  dplyr::select(qresp, alpha, theta, qnum)
md2 <- na.omit(md2)
md.aov <- lmer(qresp ~ alpha*qnum +
                 theta*qnum +
                 (1|subj) , md2)
plot(md.aov)
qqPlot(resid(md.aov))
Anova(md.aov)


### do for all levels and calculate the p-values
tal <- c()
tth <- c()
pal <- c()
pth <- c()
dal <- c()
dth <- c()
for (qn in 1:5)
{
  md2 <- filter(md.plot, qnum==qn)
  ## no 2-way int of qnum with hrv <- remove
  md.aov <- lmer(qresp ~ alpha +
                   theta + (1|subj) , md2)
  saov <- summary(md.aov)

  tal[qn] <- saov$coefficients[2,4]
  tth[qn] <- saov$coefficients[3,4]

  pal[qn] <- saov$coefficients[2,5]
  pth[qn] <- saov$coefficients[3,5]
  
  dal[qn] <- saov$coefficients[2,3]
  dth[qn] <- saov$coefficients[3,3]  
}

qname <- c("DL0","DL1","DL2","DL3","DL4")
tvalues <- data.frame(tal,tth,qname)
pvalues <- data.frame(pal,pth,qname)
colnames(tvalues) <- c("Alpha",
                       "Theta","qn")
tvalues <- tvalues %>% pivot_longer(cols = 1:2, 
                                    values_to = "tvalue",
                                    names_to = "Frequency")
pvalues <- pvalues %>% pivot_longer(cols = 1:2, 
                                    values_to = "pvalue",
                                    names_to = "Frequency")
pvalues$pvalue[pvalues$pvalue<.05] = 2
pvalues$pvalue[pvalues$pvalue<1] = 1

modtval <- ggplot(tvalues,aes(x=qn, y=tvalue, colour=Frequency, group=Frequency)) + 
  geom_line(linetype="dashed", size = .8) +
  geom_point(size = 3.5, aes(shape = factor(pvalues$pvalue)), show.legend = FALSE)+ 
  geom_point(size=2) +
  coord_cartesian(ylim = c(-5, 5)) +
  scale_y_continuous(breaks=seq(-5,5,by=2.5))+
  xlab("Depth Level") + ylab("t-value") + 
  theme_pubclean() + labs(tag = "B") +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = font_size, face = "bold"))
  # + scale_color_grey(start = .2, end = .6)


###### figure 2
grid.arrange(grobs = list(modred, modtval), nrow = 2) # 3.35 x 6 inches



#### model at each channel for topoplot
pval = array(0, c(nchans,2))
tval = array(0, c(nchans,2))
fixef= array(0, c(nchans,8))
for (nch in 1:nchans)
{
  md2 <- filter(medriv_data, 
                block=="BL" | block=="CH" | timefrend2 < lasttime,
                chans==nch) %>% 
    group_by(group, subj, block) %>%
    summarise(alpha = mean(alpha, na.rm = TRUE),
              theta = mean(theta, na.rm = TRUE),
              q1 = mean(q1, na.rm = TRUE),
              q2 = mean(q2, na.rm = TRUE),
              q3 = mean(q3, na.rm = TRUE),
              q4 = mean(q4, na.rm = TRUE),
              q5 = mean(q5, na.rm = TRUE))
  md2 <- na.omit(md2)
  
  # put all the questionnaires in a single column
  md2 <- md2 %>%
    pivot_longer(
      cols = starts_with("q"),
      names_to = "qnum",
      names_prefix = "q",
      values_to = "qresp",
      values_drop_na = TRUE
    )
  
  # md3 <- filter(md2, qnum==qn)
  md.aov <- lmer(qresp ~ theta*qnum + alpha*qnum+(1|subj), md2)
  # aov <- anova(md.aov)
  saov <- summary(md.aov)
  aov <- Anova(md.aov)
  pval[nch,] = aov$`Pr(>Chisq)`[4:5]
  tval[nch,] = aov$Chisq[4:5]
  fixef[nch,]= saov$coefficients[8:15,4]
}
write.table(c(tval,pval), row.names = FALSE, col.names = FALSE, quote = FALSE, file = 
              '~/OneDrive/Projects/Experiments/Meditation_Rivalry/Data/medpheno_alpha_theta_interaction_all.csv')


pval = array(0, c(nchans,2))
tval = array(0, c(nchans,2))
fixef= array(0, c(nchans,8))
for (nch in 1:nchans)
{
  md2 <- filter(medriv_data, 
                block=="BL" | block=="CH" | timefrend2 < lasttime,
                group=="Meditators",
                chans==nch) %>% 
    group_by(group, subj, block) %>%
    summarise(alpha = mean(alpha, na.rm = TRUE),
              theta = mean(theta, na.rm = TRUE),
              q1 = mean(q1, na.rm = TRUE),
              q2 = mean(q2, na.rm = TRUE),
              q3 = mean(q3, na.rm = TRUE),
              q4 = mean(q4, na.rm = TRUE),
              q5 = mean(q5, na.rm = TRUE))
  md2 <- na.omit(md2)
  
  # put all the questionnaires in a single column
  md2 <- md2 %>%
    pivot_longer(
      cols = starts_with("q"),
      names_to = "qnum",
      names_prefix = "q",
      values_to = "qresp",
      values_drop_na = TRUE
    )
  
  # md3 <- filter(md2, qnum==qn)
  md.aov <- lmer(qresp ~ theta*qnum + alpha*qnum+(1|subj), md2)
  saov <- summary(md.aov)
  aov <- Anova(md.aov)
  pval[nch,] = aov$`Pr(>Chisq)`[4:5]
  tval[nch,] = aov$Chisq[4:5]
  fixef[nch,]= saov$coefficients[8:15,4]
}
write.table(c(tval,pval), row.names = FALSE, col.names = FALSE, quote = FALSE, file = 
            '~/OneDrive/Projects/Experiments/Meditation_Rivalry/Data/medpheno_alpha_theta_interaction_ltm.csv')

### for all channels look at the interaction between npm and depth-level calculate the p-values
pval = array(0, c(nchans,2))
tval = array(0, c(nchans,2))
fixef= array(0, c(nchans,8))
for (nch in 1:nchans)
{
  md2 <- filter(medriv_data, 
                block=="BL" | block=="CH" | timefrend2 < lasttime,
                group=="Controls",
                chans==nch) %>% 
    group_by(group, subj, block) %>%
    summarise(alpha = mean(alpha, na.rm = TRUE),
              theta = mean(theta, na.rm = TRUE),
              q1 = mean(q1, na.rm = TRUE),
              q2 = mean(q2, na.rm = TRUE),
              q3 = mean(q3, na.rm = TRUE),
              q4 = mean(q4, na.rm = TRUE),
              q5 = mean(q5, na.rm = TRUE))
  md2 <- na.omit(md2)
  
  # put all the questionnaires in a single column
  md2 <- md2 %>%
    pivot_longer(
      cols = starts_with("q"),
      names_to = "qnum",
      names_prefix = "q",
      values_to = "qresp",
      values_drop_na = TRUE
    )
  
  # md3 <- filter(md2, qnum==qn)
  md.aov <- lmer(qresp ~ theta*qnum + alpha*qnum+(1|subj), md2)
  saov <- summary(md.aov)
  aov <- Anova(md.aov)
  pval[nch,] = aov$`Pr(>Chisq)`[4:5]
  tval[nch,] = aov$Chisq[4:5]
  fixef[nch,]= saov$coefficients[8:15,4]
}
write.table(c(tval,pval), row.names = FALSE, col.names = FALSE, quote = FALSE, file = 
            '~/OneDrive/Projects/Experiments/Meditation_Rivalry/Data/medpheno_alpha_theta_interaction_ctl.csv')



#####
## exploratory analysis for physiological measures

lasttime = 15*60

md.plot <- filter(medriv_data
                  ,chans %in% which(!(channames %in% c("FT9", "FT10", "TP9", "TP10", "FP1", "Fp2")))
                  , block=="BL" | block=="CH" | timefrend2 < lasttime 
) %>% 
  group_by(group, subj, block, subjcode) %>%
  dplyr::summarise(alpha = mean(alpha, na.rm = TRUE),
            beta = mean(beta, na.rm = TRUE),
            theta = mean(theta, na.rm = TRUE),
            gamma = mean(gamma, na.rm = TRUE),
            hrv = mean(hrv, na.rm = TRUE)/1,
            hrate = mean(hrate, na.rm = TRUE)/1,
            rrate = mean(rrate, na.rm = TRUE)/1,
            q1 = mean(q1, na.rm = TRUE),
            q2 = mean(q2, na.rm = TRUE),
            q3 = mean(q3, na.rm = TRUE),
            q4 = mean(q4, na.rm = TRUE),
            q5 = mean(q5, na.rm = TRUE))

# put all the questionnaires in a single column
md.plot <- md.plot %>%
  pivot_longer(
    cols = starts_with("q"),
    names_to = "qnum",
    names_prefix = "q",
    values_to = "qresp"
  )


md.aov <- lmer(qresp ~ 
                 hrv*qnum*group+hrate*qnum*group+rrate*qnum*group
               + (1|subj) , md.plot, REML = FALSE)
Anova(md.aov)

m2 <- update(md.aov, ~.-qnum:group:rrate)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-rrate:group)
md.aov <- m2 # accept the reduced model
Anova(md.aov)

m2 <- update(md.aov, ~.-qnum:group:hrv)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
m2 <- update(md.aov, ~.-hrv:group)
md.aov <- m2 # accept the reduced model
Anova(md.aov)

m2 <- update(md.aov, ~.-qnum:group:hrate)
anova(md.aov,m2)
# md.aov <- m2 # don't accept the reduced model

m2 <- update(md.aov, ~.-qnum:hrv)
anova(md.aov,m2)

m2 <- update(md.aov, ~.-qnum:rrate)
anova(md.aov,m2)


md.aov <- lmer(qresp ~ rrate*qnum+qnum*hrv+hrate*qnum*group
               + (1|subj) , md.plot, REML = TRUE)
plot(md.aov)
qqPlot(resid(md.aov))
Anova(md.aov)

dhrm <- c()
dhrc <- c()
thrm <- c()
thrc <- c()
phrm <- c()
phrc <- c()
for (qn in 1:5)
{
  md2 <- filter(md.plot, qnum==qn, group=="Meditators")
  md.aov <- lmer(qresp ~ hrate + (1|subj) , md2)
  # aov <- Anova(md.aov)
  saov <- summary(md.aov)
  dhrm[qn] <- saov$coefficients[2,3]
  thrm[qn] <- saov$coefficients[2,4]
  phrm[qn] <- saov$coefficients[2,5]
  
  md2 <- filter(md.plot, qnum==qn, group=="Controls")
  md.aov <- lmer(qresp ~ hrate + (1|subj) , md2)
  # aov <- Anova(md.aov)
  saov <- summary(md.aov)
  dhrc[qn] <- saov$coefficients[2,3]
  thrc[qn] <- saov$coefficients[2,4]
  phrc[qn] <- saov$coefficients[2,5]
}
phrm <- p.adjust(phrm) # fdr correction
phrc <- p.adjust(phrc) 

qname <- c("DL0","DL1","DL2","DL3","DL4")
tvalues <- data.frame(thrm,thrc,qname)
pvalues <- data.frame(phrm,phrc,qname)
colnames(tvalues) <- c("Meditators",
                       "Controls","qn")
tvalues <- tvalues %>% pivot_longer(cols = 1:2, 
                                    values_to = "tval",
                                    names_to = "HR")
p1 <- pvalues %>% pivot_longer(cols = 1:2, 
                                    values_to = "pvalue",
                                    names_to = "HR")
p1$pvalue[p1$pvalue<.05] = 2
p1$pvalue[p1$pvalue<1] = 1

hrmod <- ggplot(tvalues,aes(x=qn, y=tval, colour=HR, group=HR)) + 
  geom_line(linetype="dashed", size=.8) +
  geom_point(size = 3.5, aes(shape = factor(p1$pvalue)), show.legend = FALSE)+ 
  coord_cartesian(ylim = c(-3, 3)) + theme_pubclean() + 
  geom_point(size=2) +
  scale_y_continuous(breaks=seq(-3,3,by=1.5))+
  xlab("Depth Level") + ylab("t-value") + labs(tag = "A") +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = font_size, face = "bold")) #+ 
  # scale_color_grey(start = .2, end = .6)


drrm <- c()
dhvm <- c()
trrm <- c()
thvm <- c()
prrm <- c()
phvm <- c()
for (qn in 1:5)
{
  md2 <- filter(md.plot, qnum==qn)
  md.aov <- lmer(qresp ~ hrv + rrate + (1|subj) , md2)
  saov <- summary(md.aov)
  dhvm[qn] <- saov$coefficients[2,3]
  drrm[qn] <- saov$coefficients[3,3]

  thvm[qn] <- saov$coefficients[2,4]
  trrm[qn] <- saov$coefficients[3,4]

  phvm[qn] <- saov$coefficients[2,5]
  prrm[qn] <- saov$coefficients[3,5]
}
prrm <- p.adjust(prrm)
phvm <- p.adjust(phvm)
# 
# qname <- c("DL0","DL1","DL2","DL3","DL4")
# tvalues <- data.frame(thvm,trrm,qname)
# pvalues <- data.frame(phvm,prrm,qname)
# colnames(tvalues) <- c("HRV",  "RR", "qn")
# tvalues <- tvalues %>% pivot_longer(cols = 1:2, 
#                                     values_to = "tval",
#                                     names_to = "Physio")
# p2 <- pvalues %>% pivot_longer(cols = 1:2, 
#                                     values_to = "pvalue",
#                                     names_to = "Physio")
# p2$pvalue[p2$pvalue<.05] = 2
# p2$pvalue[p2$pvalue<=1] = 1
# 
# hrvrrmod <- ggplot(tvalues,aes(x=qn, y=tval, colour=Physio, group=Physio)) + 
#   geom_line(linetype="dashed", size=.8) +
#   geom_point(size = 3.5, aes(shape = factor(p2$pvalue)), show.legend = FALSE)+ 
#   geom_point(size=2) +
#   coord_cartesian(ylim = c(-4, 4)) + theme_pubclean() + 
#   scale_y_continuous(breaks=seq(-4,4,by=2))+
#   xlab("Depth Level") + ylab("t-value") + labs(tag = "B") +
#   theme(text = element_text(size=font_size), 
#         plot.caption = element_text(size = font_size),
#         plot.tag = element_text(size = font_size, face = "bold")) 
# # + scale_color_grey(start = .2, end = .6)

qname <- c("DL0","DL1","DL2","DL3","DL4")
tvalues <- data.frame(trrm,qname)
pvalues <- data.frame(prrm,qname)
colnames(tvalues) <- c( "RR", "qn")
tvalues <- tvalues %>% pivot_longer(cols = 1:1, 
                                    values_to = "tval",
                                    names_to = "RR")
p2 <- pvalues %>% pivot_longer(cols = 1:1, 
                               values_to = "pvalue",
                               names_to = "RR")
p2$pvalue[p2$pvalue<.05] = 2
p2$pvalue[p2$pvalue<=1] = 1

rrmod <- ggplot(tvalues,aes(x=qn, y=tval, colour=RR, group=RR)) + 
  geom_line(linetype="dashed", size=.8) +
  geom_point(size = 3.5, aes(shape = factor(p2$pvalue)), show.legend = FALSE)+ 
  geom_point(size=2) +
  coord_cartesian(ylim = c(-4, 4)) + theme_pubclean() + 
  scale_y_continuous(breaks=seq(-4,4,by=2))+
  xlab("Depth Level") + ylab("t-value") + labs(tag = "B") +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = font_size, face = "bold")) #+ 
  # scale_color_grey(start = .2, end = .6)

# qname <- c("DL0","DL1","DL2","DL3","DL4")
tvalues <- data.frame(thvm,qname)
pvalues <- data.frame(phvm,qname)
colnames(tvalues) <- c("HRV", "qn")
tvalues <- tvalues %>% pivot_longer(cols = 1:1, 
                                    values_to = "tval",
                                    names_to = "HRV")
p3 <- pvalues %>% pivot_longer(cols = 1:1, 
                               values_to = "pvalue",
                               names_to = "HRV")
p3$pvalue[p2$pvalue<.05] = 2
p3$pvalue[p2$pvalue<=1] = 1

hrvmod <- ggplot(tvalues,aes(x=qn, y=tval, colour=HRV, group=HRV)) + 
  geom_line(linetype="dashed", size=.8) +
  geom_point(size = 3.5, aes(shape = factor(p3$pvalue)), show.legend = FALSE)+ 
  geom_point(size=2) +
  coord_cartesian(ylim = c(-4, 4)) + theme_pubclean() + 
  scale_y_continuous(breaks=seq(-4,4,by=2))+
  xlab("Depth Level") + ylab("t-value") + labs(tag = "C") +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = font_size, face = "bold"))# + 
  # scale_color_grey(start = .2, end = .6)

###### figure 4
grid.arrange(grobs = list(hrmod, rrmod, hrvmod), ncol = 1)


####
## control analysis by removing the auditory stimulation part of the guided meditation to the control group

lasttime = 6*60

pval = array(0, c(nchans,2))
tval = array(0, c(nchans,2))
fixef = array(0, c(nchans,8))
for (nch in 1:nchans)
{
  md2 <- filter(medriv_data, 
                block=="BL" | block=="CH" | timefrend2 < lasttime,
                group=="Meditators",
                chans==nch) %>% 
    group_by(group, subj, block) %>%
    summarise(alpha = mean(alpha, na.rm = TRUE),
              theta = mean(theta, na.rm = TRUE),
              q1 = mean(q1, na.rm = TRUE),
              q2 = mean(q2, na.rm = TRUE),
              q3 = mean(q3, na.rm = TRUE),
              q4 = mean(q4, na.rm = TRUE),
              q5 = mean(q5, na.rm = TRUE))
  md2 <- na.omit(md2)
  
  # put all the questionnaires in a single column
  md2 <- md2 %>%
    pivot_longer(
      cols = starts_with("q"),
      names_to = "qnum",
      names_prefix = "q",
      values_to = "qresp",
      values_drop_na = TRUE
    )
  
  # md3 <- filter(md2, qnum==qn)
  md.aov <- lmer(qresp ~ theta*qnum + alpha*qnum+(1|subj), md2)
  saov <- summary(md.aov)
  aov <- Anova(md.aov)
  pval[nch,] = aov$`Pr(>Chisq)`[4:5]
  tval[nch,] = aov$Chisq[4:5]
  fixef[nch,]= saov$coefficients[8:15,3]
}
write.table(c(tval,pval), row.names = FALSE, col.names = FALSE, quote = FALSE, file = 
              '~/OneDrive/Projects/Experiments/Meditation_Rivalry/Data/medpheno_alpha_theta_interaction_ltm.csv')

### for all channels look at the interaction between npm and depth-level calculate the p-values
pval = array(0, c(nchans,2))
tval = array(0, c(nchans,2))
fixef = array(0, c(nchans,8))
for (nch in 1:nchans)
{
  md2 <- filter(medriv_data, 
                block=="BL" | block=="CH" | timefrend2 < lasttime,
                group=="Controls",
                chans==nch) %>% 
    group_by(group, subj, block) %>%
    summarise(alpha = mean(alpha, na.rm = TRUE),
              theta = mean(theta, na.rm = TRUE),
              q1 = mean(q1, na.rm = TRUE),
              q2 = mean(q2, na.rm = TRUE),
              q3 = mean(q3, na.rm = TRUE),
              q4 = mean(q4, na.rm = TRUE),
              q5 = mean(q5, na.rm = TRUE))
  md2 <- na.omit(md2)
  
  # put all the questionnaires in a single column
  md2 <- md2 %>%
    pivot_longer(
      cols = starts_with("q"),
      names_to = "qnum",
      names_prefix = "q",
      values_to = "qresp",
      values_drop_na = TRUE
    )
  
  # md3 <- filter(md2, qnum==qn)
  md.aov <- lmer(qresp ~ theta*qnum + alpha*qnum+(1|subj), md2)
  saov <- summary(md.aov)
  aov <- Anova(md.aov)
  pval[nch,] = aov$`Pr(>Chisq)`[4:5]
  tval[nch,] = aov$Chisq[4:5]  
  fixef[nch,]= saov$coefficients[8:15,3]
}
write.table(c(tval,pval), row.names = FALSE, col.names = FALSE, quote = FALSE, file = 
              '~/OneDrive/Projects/Experiments/Meditation_Rivalry/Data/medpheno_alpha_theta_interaction_ctl.csv')


####
## topographies without chanting

lasttime = 15*60

pval = array(0, c(nchans,2))
tval = array(0, c(nchans,2))
for (nch in 1:nchans)
{
  md2 <- filter(medriv_data, block!="CH",
                block=="BL" | timefrend2 < lasttime,
                group=="Meditators",
                chans==nch) %>% 
    group_by(group, subj, block) %>%
    summarise(alpha = mean(alpha, na.rm = TRUE),
              theta = mean(theta, na.rm = TRUE),
              q1 = mean(q1, na.rm = TRUE),
              q2 = mean(q2, na.rm = TRUE),
              q3 = mean(q3, na.rm = TRUE),
              q4 = mean(q4, na.rm = TRUE),
              q5 = mean(q5, na.rm = TRUE))
  md2 <- na.omit(md2)
  
  # put all the questionnaires in a single column
  md2 <- md2 %>%
    pivot_longer(
      cols = starts_with("q"),
      names_to = "qnum",
      names_prefix = "q",
      values_to = "qresp",
      values_drop_na = TRUE
    )
  
  # md3 <- filter(md2, qnum==qn)
  md.aov <- lmer(qresp ~ theta*qnum + alpha*qnum+(1|subj), md2)
  saov <- summary(md.aov)
  aov <- Anova(md.aov)
  pval[nch,] = aov$`Pr(>Chisq)`[4:5]
  tval[nch,] = aov$Chisq[4:5]  
}
write.table(c(tval,pval), row.names = FALSE, col.names = FALSE, quote = FALSE, file = 
              '~/OneDrive/Projects/Experiments/Meditation_Rivalry/Data/medpheno_alpha_theta_interaction_ltm.csv')

### for all channels look at the interaction between npm and depth-level calculate the p-values
pval = array(0, c(nchans,2))
tval = array(0, c(nchans,2))
for (nch in 1:nchans)
{
  md2 <- filter(medriv_data,  block!="CH",
                block=="BL" | timefrend2 < lasttime,
                group=="Controls",
                chans==nch) %>% 
    group_by(group, subj, block) %>%
    summarise(alpha = mean(alpha, na.rm = TRUE),
              theta = mean(theta, na.rm = TRUE),
              q1 = mean(q1, na.rm = TRUE),
              q2 = mean(q2, na.rm = TRUE),
              q3 = mean(q3, na.rm = TRUE),
              q4 = mean(q4, na.rm = TRUE),
              q5 = mean(q5, na.rm = TRUE))
  md2 <- na.omit(md2)
  
  # put all the questionnaires in a single column
  md2 <- md2 %>%
    pivot_longer(
      cols = starts_with("q"),
      names_to = "qnum",
      names_prefix = "q",
      values_to = "qresp",
      values_drop_na = TRUE
    )
  
  # md3 <- filter(md2, qnum==qn)
  md.aov <- lmer(qresp ~ theta*qnum + alpha*qnum+(1|subj), md2)
  saov <- summary(md.aov)
  aov <- Anova(md.aov)
  pval[nch,] = aov$`Pr(>Chisq)`[4:5]
  tval[nch,] = aov$Chisq[4:5]  
}
write.table(c(tval,pval), row.names = FALSE, col.names = FALSE, quote = FALSE, file = 
              '~/OneDrive/Projects/Experiments/Meditation_Rivalry/Data/medpheno_alpha_theta_interaction_ctl.csv')


####
## topographies only chanting

lasttime = 6*60

pval = array(0, c(nchans,2))
tval = array(0, c(nchans,2))
for (nch in 1:nchans)
{
  md2 <- filter(medriv_data, block=="CH",
                group=="Meditators",
                chans==nch) %>% 
    group_by(group, subj, block) %>%
    summarise(alpha = mean(alpha, na.rm = TRUE),
              theta = mean(theta, na.rm = TRUE),
              q1 = mean(q1, na.rm = TRUE),
              q2 = mean(q2, na.rm = TRUE),
              q3 = mean(q3, na.rm = TRUE),
              q4 = mean(q4, na.rm = TRUE),
              q5 = mean(q5, na.rm = TRUE))
  md2 <- na.omit(md2)
  
  # put all the questionnaires in a single column
  md2 <- md2 %>%
    pivot_longer(
      cols = starts_with("q"),
      names_to = "qnum",
      names_prefix = "q",
      values_to = "qresp",
      values_drop_na = TRUE
    )
  
  # md3 <- filter(md2, qnum==qn)
  md.aov <- lmer(qresp ~ theta*qnum + alpha*qnum+(1|subj), md2)
  saov <- summary(md.aov)
  aov <- Anova(md.aov)
  pval[nch,] = aov$`Pr(>Chisq)`[4:5]
  tval[nch,] = aov$Chisq[4:5]  
}
write.table(c(tval,pval), row.names = FALSE, col.names = FALSE, quote = FALSE, file = 
              '~/OneDrive/Projects/Experiments/Meditation_Rivalry/Data/medpheno_alpha_theta_interaction_ltm.csv')

### for all channels look at the interaction between npm and depth-level calculate the p-values
pval = array(0, c(nchans,2))
tval = array(0, c(nchans,2))
for (nch in 1:nchans)
{
  md2 <- filter(medriv_data,  block=="CH",
                group=="Controls",
                chans==nch) %>% 
    group_by(group, subj, block) %>%
    summarise(alpha = mean(alpha, na.rm = TRUE),
              theta = mean(theta, na.rm = TRUE),
              q1 = mean(q1, na.rm = TRUE),
              q2 = mean(q2, na.rm = TRUE),
              q3 = mean(q3, na.rm = TRUE),
              q4 = mean(q4, na.rm = TRUE),
              q5 = mean(q5, na.rm = TRUE))
  md2 <- na.omit(md2)
  
  # put all the questionnaires in a single column
  md2 <- md2 %>%
    pivot_longer(
      cols = starts_with("q"),
      names_to = "qnum",
      names_prefix = "q",
      values_to = "qresp",
      values_drop_na = TRUE
    )
  
  # md3 <- filter(md2, qnum==qn)
  md.aov <- lmer(qresp ~ theta*qnum + alpha*qnum+(1|subj), md2)
  saov <- summary(md.aov)
  aov <- Anova(md.aov)
  pval[nch,] = aov$`Pr(>Chisq)`[4:5]
  tval[nch,] = aov$Chisq[4:5]  
}
write.table(c(tval,pval), row.names = FALSE, col.names = FALSE, quote = FALSE, file = 
              '~/OneDrive/Projects/Experiments/Meditation_Rivalry/Data/medpheno_alpha_theta_interaction_ctl.csv')




#### exploratory analysis with medi

lasttime = 15*60

md.plot <- filter(medriv_data
                  ,chans %in% which(!(channames %in% c("FT9", "FT10", "TP9", "TP10", "FP1", "Fp2")))
                  , block=="BL" | block=="CH" | timefrend2 < lasttime 
                  # block is either baseline or chanting (~7 min blocks) or 
                  # take the last 7 minutes of the meditation blocks
) %>% 
  group_by(group, subj, block, subjcode) %>%
  summarise(alpha = mean(alpha, na.rm = TRUE),
            beta = mean(beta, na.rm = TRUE),
            theta = mean(theta, na.rm = TRUE),
            gamma = mean(gamma, na.rm = TRUE),
            hrv = mean(hrv, na.rm = TRUE)/10,
            hrate = mean(hrate, na.rm = TRUE)/10,
            rrate = mean(rrate, na.rm = TRUE)/10,
            qresp = mean(medi, na.rm = TRUE))

md.aov <- lmer(qresp ~ theta*group+alpha*group+beta*group+
                 gamma*group+
                 hrate*group+hrv*group+rrate*group
               + (1|subj) , md.plot, REML = FALSE)
Anova(md.aov)

# start by removing the 3-way ints and the associated 2-way ints with group
m2 <- update(md.aov, ~.-group:hrv)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced mode
Anova(md.aov)

m2 <- update(md.aov, ~.-group:alpha)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced mode
Anova(md.aov)

m2 <- update(md.aov, ~.-group:hrate)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
Anova(md.aov)

m2 <- update(md.aov, ~.-group:gamma)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced mode
Anova(md.aov)

m2 <- update(md.aov, ~.-group:rrate)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
Anova(md.aov)

m2 <- update(md.aov, ~.-group:beta)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced mode
Anova(md.aov)


m2 <- update(md.aov, ~.-gamma)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
Anova(md.aov)

m2 <- update(md.aov, ~.-beta)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced mode
Anova(md.aov)

md.aov <- lmer(qresp ~ theta*group+
                 hrate+hrv+rrate+alpha
               + (1|subj) , md.plot, REML = FALSE)
Anova(md.aov)

## now remove 2-way ints
m2 <- update(md.aov, ~.-hrv)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
Anova(md.aov)

m2 <- update(md.aov, ~.-hrate)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
Anova(md.aov)

m2 <- update(md.aov, ~.-rrate)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
Anova(md.aov)

m2 <- update(md.aov, ~.-theta:group)
anova(md.aov,m2)
md.aov <- m2 # accept the reduced model
Anova(md.aov)

m2 <- update(md.aov, ~.-group)
md.aov <- m2 # accept the reduced model
Anova(md.aov)

md.aov <- lmer(qresp ~ alpha + theta
               + (1|subj) , md.plot)
Anova(md.aov)
summary(md.aov)


md.aov <- lmer(theta ~ alpha
               + (1|subj) , md.plot)
Anova(md.aov)
summary(md.aov)
plot_model(md.aov, type = "pred")
