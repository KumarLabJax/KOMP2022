
setwd(file.path(workDir, outDir, out1Dir))

################### EXPORT THE DATA
		#write.csv(DataSet, paste(DataName, "Out.csv", sep = ""))

##########################################################################################
######### 	Multiple plot function
#########   http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/

			# 
			#
			# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
			# - cols:   Number of columns in layout
			# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
			#
			# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
			# then plot 1 will go in the upper left, 2 will go in the upper right, and
			# 3 will go all the way across the bottom.
			#
			multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
			  library(grid)

			  # Make a list from the ... arguments and plotlist
			  plots <- c(list(...), plotlist)

			  numPlots = length(plots)

			  # If layout is NULL, then use 'cols' to determine layout
			  if (is.null(layout)) {
				# Make the panel
				# ncol: Number of columns of plots
				# nrow: Number of rows needed, calculated from # of cols
				layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
								ncol = cols, nrow = ceiling(numPlots/cols))
			  }

			 if (numPlots==1) {
				print(plots[[1]])

			  } else {
				# Set up the page
				grid.newpage()
				pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

				# Make each plot, in the correct location
				for (i in 1:numPlots) {
				  # Get the i,j matrix positions of the regions that contain this subplot
				  matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

				  print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
												  layout.pos.col = matchidx$col))
				}
			  }
			}
######################### CONTROL OVER TIME #################################################################


	pdf(paste("Control_",DataName,"_byDate.pdf", sep = ""), 14,4)			
			for (i in PhenoList) {
				
							control$y <- control[,i] # necessary for this to work in ggplot - Petr S. 
							control$x <- control[,TestDates]
							# calculate control mean and sd for each factor
							controlMean <- mean(control$y, na.rm = T)
							controlSD <- sd(control$y, na.rm = T)
							
							control_lwrIQR_3SD <- quantile(control$y, probs=c(.25, .75), na.rm = T)[1] - 3*controlSD
							control_uprIQR_3SD <- quantile(control$y, probs=c(.25, .75), na.rm = T)[2] + 3*controlSD

				
								plot1 <- ggplot() +
								geom_point(aes(x = x,y = y,colour = Sex),data=control,shape = 16,size = 1.0) +
								geom_line(aes(x = x,y = y),data=control,size = 1.0,fun.data = mean_sdl,stat = 'summary') +
								scale_colour_brewer(guide = guide_legend(),palette = 'Paired') + labs(title = paste("ControlOnly_",i, sep = "")) +
								#theme(axis.text.x = element_text(vjust = 0.5,hjust = 1.0,angle = 90.0)) + 
								theme_bw() + #scale_x_datetime(breaks = "3 month") +
								  geom_hline(yintercept = controlMean, color = "grey", size = 1) + 
								  geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + 
								  geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1) +
								  geom_hline(yintercept = controlMean + 3*controlSD, color = "antiquewhite3", size = 1) +
								  geom_hline(yintercept = controlMean - 3*controlSD, color = "antiquewhite3", size = 1) +
								  geom_hline(yintercept = control_lwrIQR_3SD, color = "red", size = 1) +
								  geom_hline(yintercept = control_uprIQR_3SD, color = "red", size = 1) 
								
								print(plot1)

				for (j in ExperimentalFactors) {
				control$z <- control[,j]
				
				plot2 <- ggplot() +
				geom_point(aes(x = x,y = y,colour = z),data=control,shape = 16,size = 1.0) +
				geom_line(aes(x = x,y = y, colour = z),data=control,size = 1.0,fun.data = mean_sdl,stat = 'summary') +
				scale_colour_brewer(guide = guide_legend(),palette = 'Paired') + labs(title = paste("ControlOnly_",i, sep = "")) +
				#theme(axis.text.x = element_text(vjust = 0.5,hjust = 1.0,angle = 90.0)) + 
				theme_bw() + #scale_x_datetime(breaks = "3 month") +
				geom_hline(yintercept = controlMean, color = "grey", size = 1) + geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1)

				print(plot2)
				
				}
			}
	dev.off()
###########################  ALL OVER TIME  ###############################################################


	pdf(paste("ALL_",DataName,"_byDate.pdf", sep = ""), 14,4)			
		for (i in PhenoList) {
				DataSet$y <- DataSet[,i] # necessary for this to work in ggplot - Petr S. 
				DataSet$x <- DataSet[,TestDates]
				# calculate control mean and sd for each factor
				controlMean <- mean(DataSet$y, na.rm = T)
				controlSD <- sd(DataSet$y, na.rm = T)
				
				control_lwrIQR_3SD <- quantile(DataSet$y, probs=c(.25, .75), na.rm = T)[1] - 3*controlSD
				control_uprIQR_3SD <- quantile(DataSet$y, probs=c(.25, .75), na.rm = T)[2] + 3*controlSD
			
				## the angle x-label did not work!
				plot1 <- ggplot() +
				geom_point(aes(x = x,y = y,colour = Sex),data=DataSet,shape = 16,size = 1.0) +
				geom_line(aes(x = x,y = y),data=DataSet,size = 1.0,fun.data = mean_sdl,stat = 'summary') +
				scale_colour_brewer(guide = guide_legend(),palette = 'Paired') + labs(title = paste("ALL_",i, sep = "")) +
				#theme(axis.text.x = element_text(vjust = 0.5,hjust = 1.0,angle = 90.0)) + 
				theme_bw() + #scale_x_datetime(breaks = "3 month") +
				  geom_hline(yintercept = controlMean, color = "grey", size = 1) + 
				  geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + 
				  geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1) +
				  geom_hline(yintercept = controlMean + 3*controlSD, color = "antiquewhite3", size = 1) +
				  geom_hline(yintercept = controlMean - 3*controlSD, color = "antiquewhite3", size = 1) +
				  geom_hline(yintercept = control_lwrIQR_3SD, color = "red", size = 1) +
				  geom_hline(yintercept = control_uprIQR_3SD, color = "red", size = 1) +
				  guides(guide = guide_legend(title = i))
				
				
			for (j in ExperimentalFactors) {
				DataSet$z <- DataSet[,j]


				plot2 <- ggplot() +
				geom_point(aes(x = x,y = y,colour = z),data=DataSet,shape = 16,size = 1.0) +
				geom_line(aes(x = x,y = y, colour = z),data=DataSet,size = 1.0,fun.data = mean_sdl,stat = 'summary') +
				scale_colour_brewer(guide = guide_legend(),palette = 'Paired') + labs(title = paste("ALL_",i, sep = "")) +
				#theme(axis.text.x = element_text(vjust = 0.5,hjust = 1.0,angle = 90.0)) + 
				theme_bw() + #scale_x_datetime(breaks = "3 month") +
				geom_hline(yintercept = controlMean, color = "grey", size = 1) + geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1)
				print(plot2)
						}
					}
		dev.off()
		
###########################  CONTROL OVER TIME + YEAR ###############################################################
		pdf(paste("Control_",DataName,"_byYear.pdf", sep = ""), 12,10)			
			for (i in PhenoList) {
				control$y <- control[,i] # necessary for this to work in ggplot - Petr S. 
				control$x <- control[, DayOfYear]
				control$a <- as.factor(control[, Year])
				control$z <- control[,j]
				# calculate control mean and sd for each factor
				controlMean <- mean(control$y, na.rm = T)
				controlSD <- sd(control$y, na.rm = T)

				## the angle x-label did not work!
				plot1 <- ggplot() +
				geom_point(aes(x = x,y = y),data=subset(control,!is.na(a)),shape = 16,size = 1.0) +
				geom_line(aes(x = x,y = y,colour = a),data=subset(control,!is.na(a)),size = 1.0,fun.data = mean_sdl,stat = 'summary') +
				scale_colour_brewer(guide = guide_legend(),palette = 'Paired') + labs(title = paste("ControlOnly_",i, sep = "")) +
				theme(axis.text.x = element_text(vjust = 0.5,hjust = 1.0,angle = 90.0)) + 
				theme_bw() +
				geom_hline(yintercept = controlMean, color = "grey", size = 1) + geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1)
				print(plot1)

				plot1 <- plot1 + facet_grid(facets = a ~ .)
				print(plot1)
				
				}
			
	dev.off()	
###########################  ALL OVER TIME + YEAR ###############################################################
		pdf(paste("ALL_",DataName,"_byYear.pdf", sep = ""), 12,10)			
			for (i in PhenoList) {
				DataSet$y <- DataSet[,i] # necessary for this to work in ggplot - Petr S. 
				DataSet$x <- DataSet[, DayOfYear]
				DataSet$a <- as.factor(DataSet[, Year])
				# calculate control mean and sd for each factor
				controlMean <- mean(DataSet$y, na.rm = T)
				controlSD <- sd(DataSet$y, na.rm = T)

				## the angle x-label did not work!
				plot1 <- ggplot() +
				geom_point(aes(x = x,y = y),data=subset(DataSet,!is.na(a)),shape = 16,size = 1.0) +
				geom_line(aes(x = x,y = y,colour = a),data=subset(DataSet,!is.na(a)),size = 1.0,fun.data = mean_sdl,stat = 'summary') +
				scale_colour_brewer(guide = guide_legend(),palette = 'Paired') + labs(title = paste("ControlOnly_",i, sep = "")) +
				theme(axis.text.x = element_text(vjust = 0.5,hjust = 1.0,angle = 90.0)) + 
				theme_bw() +
				geom_hline(yintercept = controlMean, color = "grey", size = 1) + geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1)
				print(plot1)

				plot1 <- plot1 + facet_grid(facets = a ~ .)
				print(plot1)
				
				}
			
	dev.off()	

# ###########################  CONTROL OVER TIME of DAY ###############################################################
		# pdf(paste("Control_",DataName,"_byTimeOfDay.pdf", sep = ""), 12,4)			
			# for (i in PhenoList) {
				# control $y <- control[,i] # necessary for this to work in ggplot - Petr S. 
				# control $x <- control[, TestTime]
				# # calculate control mean and sd for each factor
				# controlMean <- mean(control$y, na.rm = T)
				# controlSD <- sd(control$y, na.rm = T)

				# ## the angle x-label did not work!
				# plot1 <- ggplot() +
				# geom_point(aes(x = x,y = y, color = Sex),data=subset(control,!is.na(x)),shape = 16,size = 1.0) +
				# stat_smooth(aes(x = x,y = y,colour = Sex, fill = Sex),data=subset(control,!is.na(x)),method = loess,formula = 'y ~ x',fullrange = TRUE,na.rm = TRUE) +
				# scale_colour_brewer(guide = guide_legend(),palette = 'Paired') + labs(title = paste("ControlOnly_",i, sep = "")) +
				# theme(axis.text.x = element_text(vjust = 0.5,hjust = 1.0,angle = 90.0)) + 
				# theme_bw() +
				# geom_hline(yintercept = controlMean, color = "grey", size = 1) + geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1)

				# print(plot1)
				
				# }
			
	# dev.off()	
# ###########################  ALL OVER TIME OF DAY ###############################################################
		# pdf(paste("ALL_",DataName,"_byTimeOfDay.pdf", sep = ""), 12,4)			
			# for (i in PhenoList) {
				# DataSet$y <- DataSet[,i] # necessary for this to work in ggplot - Petr S. 
				# DataSet$x <- DataSet[, TestTime]
				# # calculate control mean and sd for each factor
				# controlMean <- mean(DataSet$y, na.rm = T)
				# controlSD <- sd(DataSet$y, na.rm = T)

				# ## the angle x-label did not work!
				# plot1 <- ggplot() +
				# geom_point(aes(x = x,y = y, color = Sex),data=subset(DataSet,!is.na(x)),shape = 16,size = 1.0) +
				# stat_smooth(aes(x = x,y = y,colour = Sex, fill = Sex),data=subset(DataSet,!is.na(x)),method = loess,formula = 'y ~ x',fullrange = TRUE,na.rm = TRUE) +
				# scale_colour_brewer(guide = guide_legend(),palette = 'Paired') + labs(title = paste("ALL_",i, sep = "")) +
				# theme(axis.text.x = element_text(vjust = 0.5,hjust = 1.0,angle = 90.0)) + 
				# theme_bw() +
				# geom_hline(yintercept = controlMean, color = "grey", size = 1) + geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1)

				# print(plot1)
				
				# }
			
	# dev.off()	

##############	CONTROL DENSITY PLOTS	############################################################################

	pdf(paste("Control_",DataName,"_DensityPlot.pdf", sep = ""), 10,6)
			for (i in PhenoList) {
			control$y <- control[,i]

					hist0 <- ggplot() +
					theme_bw() +
					scale_fill_brewer(guide = guide_legend(),palette = 'Paired') +
					geom_histogram(aes(x = y),data=control)

					hist1 <- ggplot() +
					geom_density(aes(x = y,y = ..density..),data=control) +
					theme_bw() + labs(title = paste("ControlOnly_",i)) +
					theme(plot.title=element_text(size=8))

					hist2 <- ggplot() +
					geom_density(aes(x = y,y = ..density..,fill = Sex),data=control,alpha = 0.2784) +
					theme_bw() +
					scale_fill_brewer(guide = guide_legend(),palette = 'Paired')  + 
					theme(legend.position="none")

					hist3 <- hist2 + facet_grid(facets = Sex ~ .)

					hist4 <- qplot(sample = control$y, stat = "qq", na.rm = TRUE)

					multiplot(hist1, hist0, hist2, hist3, hist4, cols = 3)

			}
	dev.off()
##############	ALL DENSITY PLOTS	############################################################################

	pdf(paste("ALL_",DataName,"_DensityPlot.pdf", sep = ""), 10,6)
			for (i in PhenoList) {
			DataSet$y <- DataSet[,i]

					hist0 <- ggplot() +
					theme_bw() +
					scale_fill_brewer(guide = guide_legend(),palette = 'Paired') +
					geom_histogram(aes(x = y),data=DataSet)

					hist1 <- ggplot() +
					geom_density(aes(x = y,y = ..density..),data=DataSet) +
					theme_bw() + labs(title = paste("ALL_Only_",i)) +
					theme(plot.title=element_text(size=8))

					hist2 <- ggplot() +
					geom_density(aes(x = y,y = ..density..,fill = Sex),data=DataSet,alpha = 0.2784) +
					theme_bw() +
					scale_fill_brewer(guide = guide_legend(),palette = 'Paired')  + 
					theme(legend.position="none")

					hist3 <- hist2 + facet_grid(facets = Sex ~ .)

					hist4 <- qplot(sample = DataSet$y, stat = "qq")

					multiplot(hist1, hist0, hist2, hist3, hist4, cols = 3)

			}
	dev.off()
#################  ALL BOX PLOT OVER FACTORS #########################################################################

	pdf(paste("ALL_",DataName,"_BoxPlot.pdf", sep = ""), 6,4)			
			for (j in ExperimentalFactors2) {
				for (i in PhenoList) {

				DataSet$y <- DataSet[,i] # necessary for this to work in ggplot - Petr S. 
				DataSet$z <- DataSet[,j]

				# calculate control mean and sd for each factor
				controlMean <- mean(DataSet$y, na.rm = T)
				controlSD <- sd(DataSet$y, na.rm = T)

				plot1 <- ggplot() +
				geom_point(aes(x = reorder(z, y, FUN = mean, na.rm=TRUE),y = y, na.rm=TRUE),data=subset(DataSet, !is.na(x)), shape = 21,position = position_jitter(width = 0.2)) +
				geom_boxplot(aes(y = y,x = reorder(z, y, FUN = mean, na.rm=TRUE), na.rm=TRUE),data=subset(DataSet, !is.na(x)),na.rm = TRUE,outlier.size = NA, colour = "Blue") +
				theme_bw() + labs(title = i) +
				geom_hline(yintercept = controlMean, color = "grey", size = 1) + geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1)
				print (plot1)

				plot2 <- ggplot() +
				geom_point(aes(x = reorder(z, y, FUN = mean, na.rm=TRUE),y = y, na.rm=TRUE),data=subset(DataSet, !is.na(x)),shape = 21,position = position_jitter(width = 0.2)) +
				geom_crossbar(aes(y = y,x = reorder(z, y, FUN = mean, na.rm=TRUE), na.rm=TRUE),data=subset(DataSet, !is.na(x)),na.rm = TRUE, colour = "Blue", fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') +
				theme_bw() + labs(title = i) +
				geom_hline(yintercept = controlMean, color = "grey", size = 1) + geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1)
				print (plot2)
				
				plot1 <- ggplot() +
				geom_point(aes(x = z,y = y, na.rm=TRUE),data=subset(DataSet, !is.na(x)), shape = 21,position = position_jitter(width = 0.2)) +
				geom_boxplot(aes(y = y, x = z, na.rm=TRUE),data=subset(DataSet, !is.na(x)),na.rm = TRUE,outlier.size = NA, colour = "Blue") +
				theme_bw() + labs(title = i) +
				geom_hline(yintercept = controlMean, color = "grey", size = 1) + geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1)
				print (plot1)

				plot2 <- ggplot() +
				geom_point(aes(x = z,y = y, na.rm=TRUE),data=subset(DataSet, !is.na(x)),shape = 21,position = position_jitter(width = 0.2)) +
				geom_crossbar(aes(y = y,x = z, na.rm=TRUE),data=subset(DataSet, !is.na(x)),na.rm = TRUE, colour = "Blue", fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') +
				theme_bw() + labs(title = i) +
				geom_hline(yintercept = controlMean, color = "grey", size = 1) + geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1)
				print (plot2)


				}
			}
	dev.off()	

######################  CONTROL BOX PLOTS ONLY OVER FACTORS ####################################################################
	pdf(paste("Control_",DataName,"_BoxPlot.pdf", sep = ""), 6,4)			
			for (j in ExperimentalFactors2) {
				for (i in PhenoList) {

				control$y <- control[,i] # necessary for this to work in ggplot - Petr S. 
				control$z <- control[,j]

				# calculate control mean and sd for each factor
				controlMean <- mean(control$y, na.rm = T)
				controlSD <- sd(control$y, na.rm = T)

				plot1 <- ggplot() +
				geom_point(aes(x = reorder(z, y, FUN = mean, na.rm=TRUE),y = y, na.rm=TRUE),data=subset(control, !is.na(x)), shape = 21,position = position_jitter(width = 0.2)) +
				geom_boxplot(aes(y = y,x = reorder(z, y, FUN = mean, na.rm=TRUE), na.rm=TRUE),data=subset(control, !is.na(x)),na.rm = TRUE,outlier.size = NA, colour = "Blue") +
				theme_bw() + labs(title = i) +
				geom_hline(yintercept = controlMean, color = "grey", size = 1) + geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1)
				print (plot1)

				plot2 <- ggplot() +
				geom_point(aes(x = reorder(z, y, FUN = mean, na.rm=TRUE),y = y, na.rm=TRUE),data=subset(control, !is.na(x)),shape = 21,position = position_jitter(width = 0.2)) +
				geom_crossbar(aes(y = y,x = reorder(z, y, FUN = mean, na.rm=TRUE), na.rm=TRUE),data=subset(control, !is.na(x)),na.rm = TRUE, colour = "Blue",fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') +
				theme_bw() + labs(title = i) +
				geom_hline(yintercept = controlMean, color = "grey", size = 1) + geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1)
				print (plot2)
				
				plot1 <- ggplot() +
				geom_point(aes(x = z,y = y, na.rm=TRUE),data=subset(control, !is.na(x)), shape = 21,position = position_jitter(width = 0.2)) +
				geom_boxplot(aes(y = y, x = z, na.rm=TRUE),data=subset(control, !is.na(x)),na.rm = TRUE,outlier.size = NA, colour = "Blue") +
				theme_bw() + labs(title = i) +
				geom_hline(yintercept = controlMean, color = "grey", size = 1) + geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1)
				print (plot1)

				plot2 <- ggplot() +
				geom_point(aes(x = z,y = y, na.rm=TRUE),data=subset(control, !is.na(x)),shape = 21,position = position_jitter(width = 0.2)) +
				geom_crossbar(aes(y = y,x = z, na.rm=TRUE),data=subset(control, !is.na(x)),na.rm = TRUE, colour = "Blue",fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') +
				theme_bw() + labs(title = i) +
				geom_hline(yintercept = controlMean, color = "grey", size = 1) + geom_hline(yintercept = controlMean + controlSD, color = "grey", size = 1) + geom_hline(yintercept = controlMean - controlSD, color = "grey", size = 1)
				print (plot2)

				}
			}
	dev.off()
		
####################plot CONTROL relation to predicted BW ######################################################################
	pdf(paste("Control_",DataName,"_predictedBW.pdf", sep = ""), 6,4)			
				for (i in PhenoList) {
				control$y <- control[,i] # necessary for this to work in ggplot - Petr S. 
				control$x <- control[,predBW]
				plot1 <- ggplot() +
				geom_point(aes(x = x,y = y,colour = Sex),data=subset(control, !is.na(x)),shape = 21,size = 0.5,alpha = 0.5) +
				facet_grid(facets = . ~ Sex) +
				geom_smooth(aes(x = x,y = y),data= subset(control, !is.na(x)),method = lm) + 
				labs(title = i) +
				theme_bw() +
				labs(x = "Predicted BodyWeight (g)")
				print(plot1)
			}
					
				dev.off()


####################plot ALL relation to predicted BW	######################################################################
	pdf(paste("ALL_",DataName,"_predictedBW.pdf", sep = ""), 6,4)			
				for (i in PhenoList) {
				DataSet$y <- DataSet[,i] # necessary for this to work in ggplot - Petr S. 
				DataSet$x <- DataSet [,predBW]
				plot1 <- ggplot() +
				geom_point(aes(x = x,y = y,colour = Sex), data=subset(DataSet, !is.na(x)),shape = 21,size = 0.5,alpha = 0.5) +
				facet_grid(facets = . ~ Sex) +
				geom_smooth(aes(x = x,y = y), data=subset(DataSet, !is.na(x)),method = lm) +
				labs(title = i) +
				theme_bw() +
				labs(x = "Predicted BodyWeight (g)")
				print(plot1)
			}
				dev.off()


#################### CONTROL STATS ######################################################################


	### Remove ESCell from list since there is only one factor here
	ExperimentalFactors2 <- ExperimentalFactors2[ExperimentalFactors2 != ExperimentalFactors2[grep("_EScell",ExperimentalFactors2)]]



	sink(paste("Control_",DataName,"_Stats.txt", sep = ""))
			for (j in ExperimentalFactors2) {
				print("***********************************************************************")
				print(paste("**************","Control_",j,"_Stats.txt",sep = ""))
				print("***********************************************************************")

				for (i in PhenoList) {
					print(paste("****************   ",i))
					control$y <- control[,i] # necessary for this to work in ggplot - Petr S. 
					control$z <- control[,j]
					aov1 <- aov (y ~ z, data = control)
					print(summary(aov1))
					print(aov1)
					summary(control$y)
					cat('\n\n\n')
					}
				}
				
				print("***********************************************************************")
				print(paste("**************","Control_BW_",j,"_Stats.txt",sep = ""))
				print("***********************************************************************")

				### test for BW
				for (i in PhenoList) {
					print(paste("****************   ",i))
					control$y <- control[,i] # necessary for this to work in ggplot - Petr S. 
					control$z <- control[, predBW]
					aov1 <- aov (y ~ z, data = control)
					print(summary(aov1))
					print(aov1)
					summary(control$y)
					cat('\n\n\n')
					}
				
	sink()
#################### ALL STATS ######################################################################
	sink(paste("ALL_",DataName,"_Stats.txt", sep = ""))
			for (j in ExperimentalFactors2) {
				print("***********************************************************************")
				print(paste("**************","ALL_",j,"_Stats.txt",sep = ""))
				print("***********************************************************************")

				for (i in PhenoList) {
					print(paste("****************   ",i))
					DataSet$y <- DataSet[,i] # necessary for this to work in ggplot - Petr S. 
					DataSet$z <- DataSet[,j]
					aov1 <- aov (y ~ z, data = DataSet)
					print(summary(aov1))
					print(aov1)
					summary(DataSet$y)
					cat('\n\n\n')
					}
				}
					### test for BW
				print("***********************************************************************")
				print(paste("**************","ALL_BW_",j,"_Stats.txt",sep = ""))
				print("***********************************************************************")

				for (i in PhenoList) {
					print(paste("****************   ",i))
					DataSet $y <- DataSet[,i] # necessary for this to work in ggplot - Petr S. 
					DataSet $z <- DataSet[, predBW]
					aov1 <- aov (y ~ z, data = DataSet)
					print(summary(aov1))
					print(aov1)
					summary(control$y)
					cat('\n\n\n')
					}

	
	
	sink()



###################### 	PLOT PHENOTYPE DATA BY FAMILY	######################################################################
dir.create("PhenotypePlots", showWarnings = FALSE)					
setwd("PhenotypePlots")		
			
		for (i in PhenoList) {
		
		
			pdf(paste("Pheno_",i,".pdf", sep = ""), 5,20 )			

				DataSet$y <- DataSet[,i] # necessary for this to work in ggplot - Petr S. 
				DataSet$x <- DataSet[,StrainGeno]
				control$y <- control[,i] # necessary for this to work in ggplot - Petr S. 
				controlMean <- mean(control$y, na.rm = T)
				controlSD <- sd(control$y, na.rm = T)

				

				plot1 <- ggplot() +
				geom_point(aes(x = x,y = y, na.rm = TRUE),data= DataSet,shape = 20,colour = '#999999',size = 1.0) +
				geom_pointrange(aes(x = x,y = y, na.rm = TRUE),data= DataSet,colour = '#0000ff',fill = '#ffffff',size = 0.4, shape = 20,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary',na.rm=T) + 
				theme_bw(base_size = 4.0) + labs(title = i) +
				geom_hline(yintercept = controlMean, color  = "red") + 
				geom_hline(yintercept = controlMean + controlSD, color  = "red") +
				geom_hline(yintercept = controlMean - controlSD, color  = "red") + 
				coord_flip()
				print(plot1)
				
				# plot1 <- ggplot() +
				# geom_pointrange(aes(x = x,y = y, na.rm = TRUE),data= DataSet,colour = '#0000ff',fill = '#ffffff',size = 0.4, shape = 20, fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary',na.rm=T) + 
				# theme_bw(base_size = 4.0) + labs(title = i) +
				# geom_hline(yintercept = controlMean, color  = "red") + 
				# geom_hline(yintercept = controlMean + controlSD, color  = "red") +
				# geom_hline(yintercept = controlMean - controlSD, color  = "red") +
				# coord_flip()
				# print(plot1)
				
				plot1 <- ggplot() +
				geom_pointrange(aes(x = reorder(x, y, FUN = mean,na.rm=TRUE),y = y, na.rm = TRUE),data= DataSet,colour = '#0000ff',fill = '#ffffff',size = 0.4, shape = 20 ,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') + 
				theme_bw(base_size = 4.0) + labs(title = i) +
				geom_hline(yintercept = controlMean, color  = "red") + 
				geom_hline(yintercept = controlMean + controlSD, color  = "red") +
				geom_hline(yintercept = controlMean - controlSD, color  = "red") +
				geom_point(aes(x = x,y = y, na.rm = TRUE),data= DataSet,shape = 20,colour = '#999999',size =.3) +
				coord_flip()
				
				print(plot1)
											
				
				plot1 <- ggplot() +
				geom_pointrange(aes(x = reorder(x, y, FUN = mean,na.rm=TRUE),y = y, na.rm = TRUE),data= DataSet,colour = '#0000ff',fill = '#ffffff',size = 0.4, shape = 20 ,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') + 
				theme_bw(base_size = 4.0) + labs(title = i) +
				geom_hline(yintercept = controlMean, color  = "red") + 
				geom_hline(yintercept = controlMean + controlSD, color  = "red") +
				geom_hline(yintercept = controlMean - controlSD, color  = "red") +
				coord_flip()
				
				print(plot1)
				
				dev.off()

			}
			
	
####################plot all phenotypes one pdf	######################################################################

			pdf(paste("ALL_Pheno_",DataName,".pdf", sep = ""), 5,20)			

			for (i in PhenoList) {
		
				DataSet$y <- DataSet[,i] 
				DataSet$x <- DataSet[,StrainGeno]
				control$y <- control[,i] 
				controlMean <- mean(control$y, na.rm = T)
				controlSD <- sd(control$y, na.rm = T)
				
				# oldcode using pointrange - left lots of empties where there is only one animal tested
				#
				# plot1 <- ggplot() +
				# geom_pointrange(aes(x = reorder(x, y, FUN = mean,na.rm=TRUE),y = y, color = Sex, na.rm = TRUE),data=subset(DataSet, !is.na(y)),colour = '#0000ff',fill = '#ffffff',size = 0.4,shape = 20,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') + 
				# theme_bw(base_size = 4.0) + labs(title = i) +
				# geom_hline(yintercept = controlMean, color  = "red") + 
				# geom_hline(yintercept = controlMean + controlSD, color  = "red") +
				# geom_hline(yintercept = controlMean - controlSD, color  = "red") +
				# coord_flip()
				# print(plot1)
				# }
				# dev.off()
			
				plot1 <- ggplot() +
				geom_linerange(aes(x = reorder(x, y, FUN = mean,na.rm=TRUE),y = y, color = Sex, na.rm = TRUE),data=subset(DataSet, !is.na(y)),colour = '#0000ff',size = 0.4,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') + 
				geom_point(aes(x = reorder(x, y, FUN = mean,na.rm=TRUE),y = y, color = Sex, na.rm = TRUE),data=subset(DataSet, !is.na(y)),colour = '#0000ff',size = 0.8,shape = 1,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') +
				theme_bw(base_size = 4.0) + labs(title = i) +
				geom_hline(yintercept = controlMean, color  = "red") + 
				geom_hline(yintercept = controlMean + controlSD, color  = "red") +
				geom_hline(yintercept = controlMean - controlSD, color  = "red") +
				coord_flip()
				print(plot1)
				
			}
				dev.off()
