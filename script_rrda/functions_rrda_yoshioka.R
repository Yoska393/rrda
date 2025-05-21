
# Functions used for analysis of rrda simulations and applications. (Yoshioka et al.)
# These functions are used for facilitating the computation and make figures and tables.

rsbc_test <- function(Y,X,save,oneout=F){
	
	file_path <- file.path(save, "cv.RDS")
	results <-readRDS(file_path)
	num=length(results)
	b_l <- c()
	b_l.1se <-c()
	b_r <- c()
	b_r.1se <- c()
	time_cv <- c()
	
	for (i in 1:num) {
		b_l <- c(b_l, results[[i]]$cv_result$opt_min$lambda)
		b_l.1se <- c(b_l.1se, results[[i]]$cv_result$opt_lambda.1se$lambda)
		b_r <- c(b_r, results[[i]]$cv_result$opt_min$rank)
		b_r.1se <- c(b_r.1se, results[[i]]$cv_result$opt_rank.1se$rank)
		time_cv <- c(time_cv, results[[i]]$cv_result$et)
	}
	
	possible_ranks <- 1:max(b_r)
	rank_table <- table(factor(b_r, levels = possible_ranks))

	possible_ranks <- 1:max(b_r)
	rank_table <- table(factor(b_r.1se, levels = possible_ranks))

	set.seed(123)
	
	results_mspe <- list()
	
	
	pb <- txtProgressBar(min = 0, max = num, style = 3)
	for(i in 1:num){
		
		rn<-results[[i]]$rn
		
		cv_result<-results[[i]]$cv_result
		
		
		Y_train<-Y[-rn,]
		Y_test<-Y[rn,,drop=F]
		X_train<-X[-rn,]
		X_test<-X[rn,,drop=F]
		
		K<-cv_result$opt_min$rank
		L<-cv_result$opt_min$lambda
		K1<-cv_result$opt_rank.1se$rank
		L1<-cv_result$opt_lambda.1se$lambda
		K2<- min(dim(Y_train),dim(X_train))
		
		
		mspe <- c(NA,NA,NA,NA,NA,NA)
		cor <- c(NA,NA,NA,NA,NA,NA)
		ef <- c(NA,NA,NA,NA,NA,NA)
		ep <-  c(NA,NA,NA,NA,NA,NA)
		rms<- c(K,K1,K,K2,K,K2)
		lms<- c(L,L,L1,L,0,0)
		
		
		for (ms in 1: length(mspe)){
			
			tic()
			best_Bhat <- rrda.fit(Y = Y_train, X = X_train, nrank = rms[ms],lambda = lms[ms],
														scale.X = F,scale.Y = F)
			a<-toc(quiet = TRUE)
			ef[ms]<-a$toc-a$tic
			
			tic()
			best_Pred <- rrda.predict(Bhat = best_Bhat, X = X_test)[[1]][[1]][[1]]
			b<-toc(quiet = TRUE)
			ep[ms]<-b$toc-b$tic
			
			
			mspe[ms] <- 100 * sum((Y_test - best_Pred)^2) / (ncol(Y_test) * nrow(Y_test))
			if (!oneout){
				cor[ms] <- mean(diag(cor(Y_test, best_Pred)))
			}
		}
		
		
		results_mspe[[i]] <- list(
			rn = rn,
			mspe = mspe,
			cor = cor,
			time_fit = ef,
			time_pred = ep,
			L = L,
			K = K,
			L1 = L1,
			K1 = K1,
			K2 = K2
		)
		
		setTxtProgressBar(pb, i)
	}
	close(pb)
	
	file_path <- file.path(save, "test.RDS")
	dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
	saveRDS(results_mspe,file_path)
}


rsbc <- function(Y,X,num,maxrank,save,sample,Lambda=NULL,oneout=FALSE,nfold=5){
	
	if(is.null(Lambda)){
		n_lam <- 50
		Lambda <- 10^seq(0, 6, length.out = n_lam)
	}
	
	results <- list()
	
	pb <- txtProgressBar(min = 0, max = num, style = 3)
	for(i in 1:num){
		
		if (oneout){
			if (num>nrow(Y)){
				stop("num > n is not allowed")
			}
			rn <- i
		} else{
			rn<-sample(c(1:nrow(X)),sample,rep=F)
		}
		
		Y_train<-Y[-rn,]
		Y_test<-Y[rn,,drop=F]
		X_train<-X[-rn,]
		X_test<-X[rn,,drop=F]
		
		
		plan(multisession)
		tic()
		cv_result<- rrda.cv(Y = Y_train, X = X_train, lambda=Lambda,maxrank =maxrank,nfold = nfold,
												scale.X = F,scale.Y = F,verbose = F)
		a<-toc(quiet = TRUE)
		plan(sequential)
		
		K<-cv_result$opt_min$rank
		L<-cv_result$opt_min$lambda
		K1<-cv_result$opt_rank.1se$rank
		L1<-cv_result$opt_lambda.1se$lambda
		
		et<-a$toc-a$tic
		cv_result$et<-et
		results[[i]] <- list(
			rn = rn,
			cv_result = cv_result
		)
		
		setTxtProgressBar(pb, i)
	}
	close(pb)
	
	file_path <- file.path(save, "cv.RDS")
	dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
	saveRDS(results,file_path)
	
	rsbc_test(Y=Y,X=X,save=save,oneout=oneout)
}



fd <-function(col,d) {
	numeric_col <- as.numeric(col)
	formatted_col <- formatC(numeric_col , digits = d)
	return(formatted_col)
}

rd <-function(col,d) {
	numeric_col <- as.numeric(col)
	formatted_col <-  round(numeric_col, digits = d)
	formatted_col <- format(formatted_col, nsmall = d)
	return(formatted_col)
}



tp<-function(save,stars=F){
	if (stars){
		file_path <- file.path(save, "cv_stars.RDS")
	}else{
		file_path <- file.path(save, "cv.RDS")
	}
	
	results <-readRDS(file_path)
	
	num<-length(results)
	
	b_l <- c()
	b_l.1se <-c()
	b_r <- c()
	b_r.1se <- c()
	time_cv <- c()
	
	for (i in 1:num) {
		b_l <- c(b_l, results[[i]]$cv_result$opt_min$lambda)
		b_l.1se <- c(b_l.1se, results[[i]]$cv_result$opt_lambda.1se$lambda)
		b_r <- c(b_r, results[[i]]$cv_result$opt_min$rank)
		b_r.1se <- c(b_r.1se, results[[i]]$cv_result$opt_rank.1se$rank)
		time_cv <- c(time_cv, results[[i]]$cv_result$et)
	}
	print(paste("estimated maxrank",max(b_r)))

	possible_ranks <- 1:max(b_r)
	rank_table <- table(factor(b_r, levels = possible_ranks))
	barplot(rank_table, main = paste0("CV min Estimated rank"), xlab = "Rank", ylab = "Frequency")
	
	possible_ranks <- 1:max(b_r)
	rank_table <- table(factor(b_r.1se, levels = possible_ranks))
	barplot(rank_table, main = paste0("CV 1se Estimated rank"), xlab = "Rank", ylab = "Frequency")
	
	if (stars){
		file_path <- file.path(save, "test_stars.RDS")
	}else{
		file_path <- file.path(save, "test.RDS")
	}
	
	results_mspe <-readRDS(file_path)
	
	m<-c()
	for (i in 1:num) {
		m<-rbind(m,results_mspe[[i]]$mspe)
	}
	# Column names
	columns <- c("CV min", "Rank.1se", "Lambda.1se", "Full Rank", "Lambda 0", "Full Rank + Lambda 0")
	
	f=results_mspe[[1]]$K2
	# Row values for Rank and Lambda
	r <- c(mean(b_r), mean(b_r.1se), mean(b_r), f, mean(b_r), f)
	l <- c(mean(b_l), mean(b_l),  mean(b_l.1se), mean(b_l), 0, 0)
	
	rs<- c(sd(b_r),sd(b_r.1se),sd(b_r),0,sd(b_r),0) /sqrt(num)
	ls<- c(sd(b_l),sd(b_l),sd(b_l.1se),sd(b_l),0,0) /sqrt(num)
	
	Rank <- paste0(round(r, 2), " (", round(rs, 2), ")")
	Lambda <- paste0(round(l, 2), " (", round(ls, 2), ")")
	# Assuming 'm' is your matrix with PMSE values
	# PMSE Calculation
	mean_pmse <- apply(m, 2, mean)
	se_pmse <- apply(m, 2, sd) /sqrt(num)  # standard error = sd; divide by sqrt(nrow(m)) for actual SE if necessary
	
	# Format PMSE as "mean (standard error)"
	MSPE <- paste0(round(mean_pmse, 2), " (", round(se_pmse, 2), ")")
	
	# Combine Rank, Lambda, and PMSE into a table
	table_data <- data.frame(
		Rank = Rank,
		Lambda = Lambda,
		MSPE = MSPE
	)
	
	table_data <-t(table_data)
	
	# Add column names to the data frame
	colnames(table_data) <- columns
	
	# Create xtable
	xtable_output <- xtable(table_data)
	
	# Print xtable
	print(xtable_output)  # or 'latex' depending on the format you need
	
	
	m<-c()
	for (i in 1:num) {
		m<-rbind(m,results_mspe[[i]]$mspe)
	}
	
	apply(m,2,mean)
	apply(m,2,sd)
	
	rsig <- signif(r, 3)
	lsig <- signif(l, 3)
	
	default_par <- par(no.readonly = TRUE)
	
	par(mar = c(4, 4, 4, 11)) 
	

	boxplot(m, 
					main = "MSPE", 
					names = c("CV Min", "R.1se", "L.1se", "FR", "L0", "FR+L0"),
					col = c("skyblue", "lightgreen", "salmon", "lightblue", "lightpink", "lightgray"),outline = F
	)

	legend("topright", 
				 legend = c(
				 	bquote("Rank " ~ .(rsig[1]) ~ lambda ~ .(lsig[1])),
				 	bquote("Rank " ~ .(rsig[2]) ~ lambda ~ .(lsig[2])),
				 	bquote("Rank " ~ .(rsig[3]) ~ lambda ~ .(lsig[3])),
				 	bquote("Rank " ~ .(rsig[4]) ~ lambda ~ .(lsig[4])),
				 	bquote("Rank " ~ .(rsig[5]) ~ lambda ~ .(lsig[5])),
				 	bquote("Rank " ~ .(rsig[6]) ~ lambda ~ .(lsig[6]))
				 ),
				 fill = c("skyblue", "lightgreen", "salmon", "lightblue", "lightpink", "lightgray"), 
				 title = "Methods", 
				 xpd = TRUE,          
				 inset = c(-0.45, 0), 
				 cex = 0.9)    
	
	
	m<-c()
	for (i in 1:num) {
		m<-rbind(m,results_mspe[[i]]$cor)
	}
	
	apply(m,2,mean)
	apply(m,2,sd)
	
	rsig <- signif(r, 3)
	lsig <- signif(l, 3)
	
	par(mar = c(4, 4, 4, 11))  
	
	
	boxplot(m, 
					main = "Correlation Coefficient", 
					names = c("CV Min", "R.1se", "L.1se", "FR", "L0", "FR+L0"),
					col = c("skyblue", "lightgreen", "salmon", "lightblue", "lightpink", "lightgray"),outline = F
	)
	
	
	legend("topright", 
				 legend = c(
				 	bquote("Rank " ~ .(rsig[1]) ~ lambda ~ .(lsig[1])),
				 	bquote("Rank " ~ .(rsig[2]) ~ lambda ~ .(lsig[2])),
				 	bquote("Rank " ~ .(rsig[3]) ~ lambda ~ .(lsig[3])),
				 	bquote("Rank " ~ .(rsig[4]) ~ lambda ~ .(lsig[4])),
				 	bquote("Rank " ~ .(rsig[5]) ~ lambda ~ .(lsig[5])),
				 	bquote("Rank " ~ .(rsig[6]) ~ lambda ~ .(lsig[6]))
				 ),
				 fill = c("skyblue", "lightgreen", "salmon", "lightblue", "lightpink", "lightgray"), 
				 title = "Methods", 
				 xpd = TRUE,  
				 inset = c(-0.45, 0), 
				 cex = 0.9) 
	
	par(default_par)
}


criteria_test<- function(Y,X,num,Lambda=NULL,rn_list=NULL,lm=NULL,dir_path,maxrank=20,nfold=5,
												 criteria= c("AIC", "BIC", "GIC", "BICP", "GCV","rrda","StARS"),oneout=FALSE){
	
	# Initialize a data frame to store results for each method
	results <- data.frame(Method = character(),
												Ranks_Vectors = I(list()),
												Lambdas_Vectors = I(list()),
												Average_Rank = numeric(),
												Std_Err_Rank = numeric(),
												Average_Lambda = numeric(),
												Std_Err_Lambda = numeric(),
												Average_Computation_Time = numeric(),
												Std_Err_Computation_Time = numeric(),
												MSPE_Vectors = I(list()),
												Average_MSPE = numeric(),
												Std_Err_MSPE = numeric(),
												stringsAsFactors = FALSE)
	cv_res <- list()
	# Loop over each model type
	for (mot in criteria) {
		#for (mot in c("AIC", "BIC","GIC", "BICP", "GCV")) {
		
		av <- c()  # Vector to store ranks for this method
		al <- c()  # Vector to store optimized lambdas for this method
		mv <- c()  # Vector to store MSPE values
		computation_times <- numeric(num)
		pb <- txtProgressBar(min = 0, max = num, style = 3)
		
		for (ite in 1:num) {
			
			if (oneout){
				if (num>nrow(Y)){
					stop("num > n is not allowed")
				}
				rn <- ite
			} else{
				rn <-  rn_list[[ite]]
			}
			
			Y_train <- Y[-rn,,drop=F]
			Y_test <- Y[rn,,drop=F]
			X_train <- X[-rn,,drop=F]
			X_test <- X[rn,,drop=F]
			
			if (is.null(lm)){
				if(is.null(Lambda)){
					Lambda <- exp(seq(8, 0, length.out = 50))
				}else{
					Lambda <- Lambda
				}
			} else{
				Lambda <- lm[[ite]]
			}
			
			
			if (mot == "rrda"){
				start_time <- proc.time()
				
				
				cv_result <- rrda.cv(Y = Y_train, X = X_train,lambda = Lambda,nfold=nfold,maxrank=maxrank,
														 scale.X = F,scale.Y = F,verbose = F)
				
				end_time <- proc.time()
				et <- (end_time - start_time)[3]
				cv_result$et<-et
				cv_res[[ite]]<- list(
					rn = rn,
					cv_result = cv_result
				)
				
				optimized_lambda <- cv_result$opt_min$lambda
				optimized_rank <- cv_result$opt_min$rank
			} else if (mot == "StARS") {
				start_time <- proc.time()
				rrr_result <- rrr_ann_stars2(response = Y_train, X = X_train, lambda = Lambda)
				
				end_time <- proc.time()
				
				optimized_lambda <- rrr_result$lambda
				optimized_rank <- rrr_result$rank
			} else {
				start_time <- proc.time()
				rrr_result <- rrr2(Y = Y_train, X = X_train, penaltySVD = c("ann"), ic.type = mot, modstr =  rrr.modstr(lambda=Lambda))
				
				end_time <- proc.time()
				
				optimized_lambda <- rrr_result$lambda[which.min(rrr_result$ic)]
				optimized_rank <- rrr_result$rank
			}
			
			computation_times[ite] <- (end_time - start_time)[3]
			# Append results
			av <- c(av, optimized_rank)
			al <- c(al, optimized_lambda)
			
			if (optimized_rank < 1) {
				meansYtr <- colMeans(Y_train)
				best_Pred <- matrix(meansYtr,nrow(Y_test),ncol(Y_test), byrow = TRUE)
			} else {
				model_rank <-optimized_rank
				best_Bhat <- rrda.fit(Y = Y_train, X = X_train, nrank = model_rank, lambda = optimized_lambda,scale.X = F,scale.Y = F)
				best_Pred <- rrda.predict(Bhat = best_Bhat, X = X_test)[[1]][[1]][[1]]
			}
			

			mspe <- 100 * sum((Y_test - best_Pred)^2) / (ncol(Y_test) * nrow(Y_test))
			mv <- c(mv, mspe)  # Store MSPE
			
			# Update the progress bar
			setTxtProgressBar(pb, ite)
		}
		
		close(pb)
		
		# Calculate mean and standard error for ranks, lambdas, and MSPE
		mean_rank <- mean(av)
		std_err_rank <- sd(av) / sqrt(length(av))
		
		mean_lambda <- mean(al)
		std_err_lambda <- sd(al) / sqrt(length(al))
		
		mean_mspe <- mean(mv)
		std_err_mspe <- sd(mv) / sqrt(length(mv))
		
		# Calculate mean computation time and its standard error
		mean_computation_time <- mean(computation_times)
		std_err_computation_time <- sd(computation_times) / sqrt(length(computation_times))
		
		# Populate the results data frame
		results <- rbind(results, data.frame(Method = mot,
																				 Ranks_Vectors = I(list(av)),
																				 Lambdas_Vectors = I(list(al)),
																				 Average_Rank = mean_rank,
																				 Std_Err_Rank = std_err_rank,
																				 Average_Lambda = mean_lambda,
																				 Std_Err_Lambda = std_err_lambda,
																				 Average_Computation_Time = mean_computation_time,
																				 Std_Err_Computation_Time = std_err_computation_time,
																				 MSPE_Vectors = I(list(mv)),  # Store vector of MSPE
																				 Average_MSPE = mean_mspe,     # Store mean MSPE
																				 Std_Err_MSPE = std_err_mspe    # Store standard error of MSPE
		))
		
	}
	
	
	if (!dir.exists(dir_path)) {
		dir.create(dir_path, recursive = TRUE)
	}
	
	saveRDS(cv_res, file.path(dir_path,"cv.RDS"))
	saveRDS(results, file.path(dir_path,"criteria.RDS"))
	
	return(results)
}


criteria_fit<- function(Y,X,num,Lambda=NULL,rn_list=NULL,lm=NULL,dir_path,maxrank=20,nfold=5,
												 criteria= c("AIC", "BIC", "GIC", "BICP", "GCV","rrda","StARS"),oneout=FALSE){
	
	res_criteria<-readRDS(file.path(dir_path,"criteria.RDS"))
	# Initialize a data frame to store results for each method

	# Loop over each model type
	for (mot in criteria) {
		#for (mot in c("AIC", "BIC","GIC", "BICP", "GCV")) {
		
		av <- res_criteria[res_criteria$Method==mot,]$Ranks_Vectors [[1]]
		al <-  res_criteria[res_criteria$Method==mot,]$Lambdas_Vectors [[1]]
		mv <- c()  # Vector to store MSPE values
		computation_times <- numeric(num)
		pb <- txtProgressBar(min = 0, max = num, style = 3)
		
		for (ite in 1:num) {
			
			if (oneout){
				if (num>nrow(Y)){
					stop("num > n is not allowed")
				}
				rn <- ite
			} else{
				rn <-  rn_list[[ite]]
			}
			
			Y_train <- Y[-rn,,drop=F]
			Y_test <- Y[rn,,drop=F]
			X_train <- X[-rn,,drop=F]
			X_test <- X[rn,,drop=F]

			# Append results
			
			optimized_rank   <- av[ite]
			optimized_lambda <- al[ite]
			
			if (optimized_rank < 1) {
				meansYtr <- colMeans(Y_train)
				best_Pred <- matrix(meansYtr,nrow(Y_test),ncol(Y_test), byrow = TRUE)
			} else {
				model_rank <-optimized_rank
				best_Bhat <- rrda.fit(Y = Y_train, X = X_train, nrank = model_rank, lambda = optimized_lambda,scale.X = F,scale.Y = F)
				best_Pred <- rrda.predict(Bhat = best_Bhat, X = X_test)[[1]][[1]][[1]]
			}
			
			mspe <- 100 * sum((Y_test - best_Pred)^2) / (ncol(Y_test) * nrow(Y_test))
			mv <- c(mv, mspe)  # Store MSPE
			
			# Update the progress bar
			setTxtProgressBar(pb, ite)
		}
		
		close(pb)
		
		# Calculate mean and standard error for ranks, lambdas, and MSPE
		mean_rank <- mean(av)
		std_err_rank <- sd(av) / sqrt(length(av))
		
		mean_lambda <- mean(al)
		std_err_lambda <- sd(al) / sqrt(length(al))
		
		mean_mspe <- mean(mv)
		std_err_mspe <- sd(mv) / sqrt(length(mv))

		res_criteria[res_criteria$Method==mot,]$MSPE_Vectors <- I(list(mv))
		res_criteria[res_criteria$Method==mot,]$Average_MSPE <- mean_mspe
		res_criteria[res_criteria$Method==mot,]$Std_Err_MSPE <-std_err_mspe
			
		# Populate the results data frame
		
	}
	
	
	if (!dir.exists(dir_path)) {
		dir.create(dir_path, recursive = TRUE)
	}

	saveRDS(res_criteria, file.path(dir_path,"criteria2.RDS"))
	
	return(res_criteria)
}


criteria_b<- function(Y,X,B,num,Lambda=NULL,rn_list=NULL,lm=NULL,dir_path,maxrank=20,nfold=5,
												criteria= c("AIC", "BIC", "GIC", "BICP", "GCV","rrda","StARS"),oneout=FALSE){
	
	res_criteria<-readRDS(file.path(dir_path,"criteria.RDS"))
	# Initialize a data frame to store results for each method
	
	# Loop over each model type
	for (mot in criteria) {
		#for (mot in c("AIC", "BIC","GIC", "BICP", "GCV")) {
		
		av <- res_criteria[res_criteria$Method==mot,]$Ranks_Vectors [[1]]
		al <-  res_criteria[res_criteria$Method==mot,]$Lambdas_Vectors [[1]]
		mv <- c()  # Vector to store MSPE values
		computation_times <- numeric(num)
		pb <- txtProgressBar(min = 0, max = num, style = 3)
		
		for (ite in 1:num) {
			
			if (oneout){
				if (num>nrow(Y)){
					stop("num > n is not allowed")
				}
				rn <- ite
			} else{
				rn <-  rn_list[[ite]]
			}
			
			Y_train <- Y[-rn,,drop=F]
			Y_test <- Y[rn,,drop=F]
			X_train <- X[-rn,,drop=F]
			X_test <- X[rn,,drop=F]

			
			# Append results
			
			optimized_rank   <- av[ite]
			optimized_lambda <- al[ite]
			
			if (optimized_rank < 1) {
				model_rank <-optimized_rank
				meansYtr <- colMeans(Y_train)
				best_Bhat <- rrda.fit(Y = Y_train, X = X_train, nrank = 1, lambda = optimized_lambda,scale.X = F,scale.Y = F)
				best_Bhat<-rrda.coef(best_Bhat)[[1]][[1]]
				#best_Pred <- matrix(meansYtr,nrow(Y_test),ncol(Y_test), byrow = TRUE)
			} else {
				model_rank <-optimized_rank
				best_Bhat <- rrda.fit(Y = Y_train, X = X_train, nrank = model_rank, lambda = optimized_lambda,scale.X = F,scale.Y = F)
				best_Bhat<-rrda.coef(best_Bhat)[[1]][[1]]
				#best_Pred <- rrda.predict(Bhat = best_Bhat, X = X_test)[[1]][[1]][[1]]
			}
			
			#mspe <- 100 * sum((Y_test - best_Pred)^2) / (ncol(Y_test) * nrow(Y_test))
			#mspe <- 100 * sum((B - best_Bhat)^2) / (ncol(B) * nrow(B))
			norm<-norm(B- best_Bhat, type = "F")
			#dim(best_Bhat)
			#dim(B)
			#mspe<-cor(as.vector(best_Bhat), as.vector(B))
			mv <- c(mv, norm)  # Store MSPE
			
			# Update the progress bar
			setTxtProgressBar(pb, ite)
		}
		
		close(pb)
		
		# Calculate mean and standard error for ranks, lambdas, and MSPE
		mean_rank <- mean(av)
		std_err_rank <- sd(av) / sqrt(length(av))
		
		mean_lambda <- mean(al)
		std_err_lambda <- sd(al) / sqrt(length(al))
		
		mean_norm <- mean(mv)
		std_err_norm <- sd(mv) / sqrt(length(mv))
		
		res_criteria[res_criteria$Method==mot,]$MSPE_Vectors <- I(list(mv))
		res_criteria[res_criteria$Method==mot,]$Average_MSPE <- mean_norm
		res_criteria[res_criteria$Method==mot,]$Std_Err_MSPE <-std_err_norm
		
		colnames(res_criteria)[colnames(res_criteria) == "MSPE_Vectors"] <- "NORM_Vectors"
		colnames(res_criteria)[colnames(res_criteria) == "Average_MSPE"] <- "Average_NORM"
		colnames(res_criteria)[colnames(res_criteria) == "Std_Err_MSPE"] <- "Std_Err_NORM"
		
		
		# Populate the results data frame
		
	}
	
	
	if (!dir.exists(dir_path)) {
		dir.create(dir_path, recursive = TRUE)
	}
	
	saveRDS(res_criteria, file.path(dir_path,"criteria_b.RDS"))
	
	return(res_criteria)
}

criteria_test_all<- function(Y,X,Lambda=NULL,dir_path,maxrank=20,nfold=5,
														 criteria= c("AIC", "BIC", "GIC", "BICP", "GCV","rrda","StARS")){
	
	# Initialize a data frame to store results for each method
	results <- data.frame(Method = character(),
												#Ranks_Vectors = I(list()),
												#Lambdas_Vectors = I(list()),
												Rank = numeric(),
												#Std_Err_Rank = numeric(),
												Lambda = numeric(),
												#Std_Err_Lambda = numeric(),
												Computation_Time = numeric(),
												#Std_Err_Computation_Time = numeric(),
												#MSPE_Vectors = I(list()),
												#Average_MSPE = numeric(),
												#Std_Err_MSPE = numeric(),
												stringsAsFactors = FALSE)
	
	pb <- txtProgressBar(min = 0, max = length(criteria), style = 3)
	# Loop over each model type
	mn=0
	for (mot in criteria) {
		#for (mot in c("AIC", "BIC","GIC", "BICP", "GCV")) {
		
		av <- c()  # Vector to store ranks for this method
		al <- c()  # Vector to store optimized lambdas for this method
		#mv <- c()  # Vector to store MSPE values
		computation_times <- numeric(1)
		
		if(is.null(Lambda)){
			Lambda <- exp(seq(8, 0, length.out = 50))
		}else{
			Lambda <- Lambda
		}
		
		
		if (mot == "rrda"){
			start_time <- proc.time()
			
			
			cv_result <- rrda.cv(Y = Y, X = X,lambda = Lambda,nfold=nfold,maxrank=maxrank,
													 scale.X = F,scale.Y = F,verbose = F)
			
			end_time <- proc.time()
			
			optimized_lambda <- cv_result$opt_min$lambda
			optimized_rank <- cv_result$opt_min$rank
		} else if (mot == "StARS") {
			start_time <- proc.time()
			rrr_result <- rrr_ann_stars2(response = Y, X = X, lambda = Lambda)
			
			end_time <- proc.time()
			
			optimized_lambda <- rrr_result$lambda
			optimized_rank <- rrr_result$rank
		} else {
			start_time <- proc.time()
			rrr_result <- rrr2(Y = Y, X = X, penaltySVD = c("ann"), ic.type = mot, modstr =  rrr.modstr(lambda=Lambda))
			
			end_time <- proc.time()
			
			optimized_lambda <- rrr_result$lambda[which.min(rrr_result$ic)]
			optimized_rank <- rrr_result$rank
		}
		
		computation_times <- (end_time - start_time)[3]
		
		# Append results
		av <- c(av, optimized_rank)
		al <- c(al, optimized_lambda)
		
		
		
		mean_rank <- mean(av)
		mean_lambda <- mean(al)
		mean_computation_time <- mean(computation_times)
		
		# Populate the results data frame
		results <- rbind(results, data.frame(Method = mot,
																				 #Ranks_Vectors = I(list(av)),
																				 #Lambdas_Vectors = I(list(al)),
																				 Rank = mean_rank,
																				 #Std_Err_Rank = std_err_rank,
																				 Lambda = mean_lambda,
																				 #Std_Err_Lambda = std_err_lambda,
																				 Computation_Time = mean_computation_time
																				 #Std_Err_Computation_Time = std_err_computation_time,
																				 #MSPE_Vectors = I(list(mv)),  # Store vector of MSPE
																				 #Average_MSPE = mean_mspe,     # Store mean MSPE
																				 #Std_Err_MSPE = std_err_mspe    # Store standard error of MSPE
		))
		mn<-mn+1
		setTxtProgressBar(pb, mn)
	}
	
	close(pb)
	
	
	if (!dir.exists(dir_path)) {
		dir.create(dir_path, recursive = TRUE)
	}
	
	saveRDS(results, file.path(dir_path,"criteria_alldata.RDS"))
	
	return(results)
}



rrr_lambda <- function (Y, X, penaltySVD = c("ann"), 
												ic.type = c("GIC", "AIC", "BIC", "BICP", "GCV"), 
												df.type = c("exact", "naive"), 
												maxrank = min(dim(Y), dim(X)), 
												modstr = list(), 
												control = list()) 
{
	
	penaltySVD <- match.arg(penaltySVD)
	ic.type <- match.arg(ic.type)
	df.type <- match.arg(df.type)
	q <- ncol(Y)
	n <- nrow(Y)
	p <- ncol(X)
	if (n != nrow(X)) 
		stop("'nrow(X)' has to be equal to 'nrow(Y)'.")
	control <- do.call("rrr.control", control)
	modstr <- do.call("rrr.modstr", modstr)
	qrX <- qr(X, tol = control$qr.tol)
	C_ls <- qr.coef(qrX, Y)
	C_ls <- ifelse(is.na(C_ls), 0, C_ls)
	rX <- qrX$rank
	nrank <- min(q, rX, maxrank)
	rmax <- min(rX, q)
	XC <- qr.fitted(qrX, Y)
	svdXC <- svd(XC, nu = rmax, nv = rmax)
	A <- svdXC$v[, 1:rmax]
	Ad <- (svdXC$d[1:rmax])^2
	Ad <- Ad[Ad > control$sv.tol]
	rmax <- length(Ad)
	if (identical(penaltySVD, "ann")) {
		gamma <- modstr$gamma
		nlambda <- modstr$nlambda
		lambda <- modstr$lambda
		Adx <- Ad^((1 + gamma) / 2)
		if (is.null(lambda)) {
			lambda <- exp(seq(log(max(Adx)), log(Adx[min(nrank + 
																									 	1, rmax)]), length = nlambda))
		}
	}
	
	return(lambda)
}



get_lambda_stars <- function(dir_path, rn_list, Y, X) {
	# Initialize the 100x100 matrix to store lambda values
	lm <- list()
	
	# Loop through each index in rn_list
	for (ite in 1:100) {
		rn <- rn_list[[ite]]
		
		# Create training and testing sets
		Y_train <- Y[-rn,,drop=F]
		Y_test <- Y[rn,,drop=F]
		X_train <- X[-rn,,drop=F]
		X_test <- X[rn,,drop=F]
		
		# Run the rrr2 function with the specified penalty
		Lambda<- rrr_lambda(Y=Y_train, X=X_train, penaltySVD = c("ann"))
		
		lm[[ite]] <- Lambda
	}
	
	# Create the directory if it doesn't exist
	
	if (!dir.exists(dir_path)) {
		dir.create(dir_path, recursive = TRUE)
	}
	
	# Save the lm matrix as an RDS file
	saveRDS(lm, file.path(dir_path, "aic_lambda.RDS"))
	
	# Return the lm matrix for reference
	return(lm)
}


measure_rrda_times <- function(Y, X, dir_path) {
	# Initialize lists to store times
	cv_times <- numeric(10)
	fit_times <- numeric(10)
	predict_times <- numeric(10)
	
	# Loop for 10 iterations
	for (i in 1:10) {
		# Measure rrda.cv computation time
		plan(multisession)
		cv_times[i] <- system.time({
			cv_result <- rrda.cv(Y = Y, X = X, nfold = 5, scale.X = FALSE, scale.Y = FALSE)
		})["elapsed"]
		plan(sequential)
		
		# Extract optimal rank and lambda from cv_result
		K <- cv_result$opt_min$rank
		L <- cv_result$opt_min$lambda
		
		# Measure rrda.fit computation time
		fit_times[i] <- system.time({
			best_Bhat <- rrda.fit(Y = Y, X = X, nrank = K, lambda = L, scale.X = FALSE, scale.Y = FALSE)
		})["elapsed"]
		
		# Measure rrda.predict computation time
		predict_times[i] <- system.time({
			best_Pred <- rrda.predict(Bhat = best_Bhat, X = X)[[1]][[1]][[1]]
		})["elapsed"]
	}
	
	# Print average computation times
	cat("Average rrda.cv time:", mean(cv_times), "seconds\n")
	cat("Average rrda.fit time:", mean(fit_times), "seconds\n")
	cat("Average rrda.predict time:", mean(predict_times), "seconds\n")
	
	# Store time results in a list
	time_results <- list(
		cv_times = cv_times,
		fit_times = fit_times,
		predict_times = predict_times
	)
	
	# Save the list of times to the specified directory
	saveRDS(time_results, file.path(dir_path, "comp.RDS"))
	
	# Return the time results list for reference
	return(time_results)
}

# Function to calculate mean and standard error
mean_se <- function(x) {
	mean_x <- mean(x)
	se_x <- sd(x) / sqrt(length(x))
	c(Mean = mean_x, SE = se_x)
}

compute_time_stats <- function(dir_path) {
	# Load the saved computation times
	time_results <- readRDS(file.path(dir_path, "comp.RDS"))
	
	# Apply mean and SE calculation to each time component
	cv_stats <- mean_se(time_results$cv_times)
	fit_stats <- mean_se(time_results$fit_times)
	predict_stats <- mean_se(time_results$predict_times)
	
	# Combine the results into a data frame
	results_df <- data.frame(
		Function = c("rrda.cv", "rrda.fit", "rrda.predict"),
		Mean = c(cv_stats[1], fit_stats[1], predict_stats[1]),
		SE = c(cv_stats[2], fit_stats[2], predict_stats[2])
	)
	
	# Print the table using xtable
	print(xtable(results_df, caption = "Mean and Standard Error of Computation Times", digits = 4))
	
	# Return the results data frame for reference
	return(results_df)
}



compare<-function(save,save2,chr,num,model=1){
	
	file_path <- file.path(save, "test.RDS")
	results_mspe <-readRDS(file_path)
	
	m<-c()
	for (i in 1:num) {
		m<-rbind(m,results_mspe[[i]]$mspe)
	}
	ch_all <-m
	
	file_path <- file.path(save2, "test.RDS")
	results_mspe <-readRDS(file_path)
	
	m<-c()
	for (i in 1:num) {
		m<-rbind(m,results_mspe[[i]]$mspe)
	}
	ch <-m
	
	maxv<- max(ch[,model],ch_all[,model])
	minv<- min(ch[,model],ch_all[,model])
	plot(ch[,model],ch_all[,model],main="MSPE",
			 ylab = paste0("All to ch",chr),
			 xlab = paste0("ch",chr," to ch",chr),
			 xlim=c(minv,maxv),ylim=c(minv,maxv)
	)
	abline(0,1,col="red")
	
	ch_mspe<-ch[,model]
	ch_all_mspe<-ch_all[,model]
	
	file_path <- file.path(save, "test.RDS")
	results_mspe <-readRDS(file_path)
	
	m<-c()
	for (i in 1:num) {
		m<-rbind(m,results_mspe[[i]]$cor)
	}
	ch_all <-m
	
	file_path <- file.path(save2, "test.RDS")
	results_mspe <-readRDS(file_path)
	
	m<-c()
	for (i in 1:num) {
		m<-rbind(m,results_mspe[[i]]$cor)
	}
	ch <-m
	
	maxv<- max(ch[,model],ch_all[,model])
	minv<- min(ch[,model],ch_all[,model])
	plot(ch[,model],ch_all[,model],main="COR",ylab = paste0("All to ch",chr),xlab = paste0("ch",chr," to ch",chr),
			 xlim=c(minv,maxv),ylim=c(minv,maxv)
	)
	abline(0,1,col="red")

	res<-list(ch=ch[,model],ch_all=ch_all[,model],ch_mspe=ch_mspe,ch_all_mspe=ch_all_mspe)
	return(res)
}

# Define a function to apply for each set of parameters 1
run_simulation <- function(n, p, q, k,rdasim) {
	# Generate simulated data
	if (rdasim=="sim1"){
		simdata <- rdasim1(n = n, p = p, q = q, k = k)
	} else if (rdasim=="sim2"){
		simdata <- rdasim2(n = n, p = p, q = q, k = k)
	} else{
		stop("simulation model must be selected properly")
	}
	X <- simdata$X
	Y <- simdata$Y
	
	# Run cross-validation using rrda.cv
	cv_result <- rrda.cv(Y = Y, X = X, center.X = TRUE, center.Y = TRUE,scale.X = F,scale.Y = F,verbose = F,nfold = 5)
	rrda.plot(cv_result)
	# Extract the optimal rank (esR)
	esR <- cv_result$opt_min$rank
	
	return(esR)
}

# Define a function to apply for each set of parameters 2
rank_evaluate<- function(n_iterations,rdasim,size){
	
	# Generate randomized lists for n, p, q, and k
	
	if (size == "large"){
		n_list <- sample(50:100, n_iterations, replace = TRUE)
		p_list <- sample(100:1000, n_iterations, replace = TRUE)
		q_list <- sample(100:1000, n_iterations, replace = TRUE)
		k_list <- sample(1:10, n_iterations, replace = TRUE)
	} else if (size == "small"){
		n_list <- sample(100:200, n_iterations, replace = TRUE)
		p_list <- sample(20:100, n_iterations, replace = TRUE)
		q_list <- sample(20:100, n_iterations, replace = TRUE)
		k_list <- sample(1:10, n_iterations, replace = TRUE)
	} else if (size == "mixed"){
		n_list <- sample(c(50:200), n_iterations, replace = TRUE)
		p_list <- sample(c(20:1000), n_iterations, replace = TRUE)
		q_list <- sample(c(20:1000), n_iterations, replace = TRUE)
		k_list <- sample(1:10, n_iterations, replace = TRUE)
	} else{
		stop ("size error")
	}
	
	# Apply the simulation function to each set of n, p, q, k values using map functions
	print(rdasim)
	
	# Initialize the progress bar
	pb <- progress_bar$new(
		total = length(n_list),
		format = "  Progress [:bar] :percent in :elapsed"
	)
	
	# Define a custom function to run the simulation and update the progress bar
	run_simulation_with_progress <- function(n, p, q, k, rdasim) {
		pb$tick()  # Update the progress bar
		run_simulation(n, p, q, k, rdasim)  # Call your main simulation function
	}
	
	# Use pmap_dbl with the custom function
	esR_vector <- pmap_dbl(list(n_list, p_list, q_list, k_list, rdasim = rdasim), 
												 run_simulation_with_progress)
	
	#esR_vector <- pmap_dbl(list(n_list, p_list, q_list, k_list, rdasim = rdasim), run_simulation)
	
	data_list <- list(n_list = n_list, p_list = p_list, q_list = q_list, k_list = k_list, esR_vector = esR_vector)
	
	return(data_list)
}

rank_evaluate_50<- function(n_iterations,rdasim,size){
	
	# Generate randomized lists for n, p, q, and k
	
	if (size == "small"){
		combo <- expand.grid(
			n = c(100, 200),
			pq = c(50, 100),
			k = c(2, 5, 10)
		)
		
		# Repeat each row 50 times
		combo <- combo[rep(1:nrow(combo), each = n_iterations), ]
		
		# Extract into separate lists or vectors
		n_list <- combo$n
		p_list <- combo$pq
		q_list <- combo$pq
		k_list <- combo$k
	} else if (size == "large"){
		combo <- expand.grid(
			n = c(100, 200),
			pq = c(500, 1000),
			k = c(2, 5, 10)
		)
		# Repeat each row 50 times
		combo <- combo[rep(1:nrow(combo), each = n_iterations), ]
		
		# Extract into separate lists or vectors
		n_list <- combo$n
		p_list <- combo$pq
		q_list <- combo$pq
		k_list <- combo$k
	} else{
		stop ("size error")
	}
	
	# Apply the simulation function to each set of n, p, q, k values using map functions
	print(rdasim)
	
	# Initialize the progress bar
	pb <- progress_bar$new(
		total = length(n_list),
		format = "  Progress [:bar] :percent in :elapsed"
	)
	
	# Define a custom function to run the simulation and update the progress bar
	run_simulation_with_progress <- function(n, p, q, k, rdasim) {
		pb$tick()  # Update the progress bar
		run_simulation(n, p, q, k, rdasim)  # Call your main simulation function
	}
	
	# Use pmap_dbl with the custom function
	esR_vector <- pmap_dbl(list(n_list, p_list, q_list, k_list, rdasim = rdasim), 
												 run_simulation_with_progress)
	
	#esR_vector <- pmap_dbl(list(n_list, p_list, q_list, k_list, rdasim = rdasim), run_simulation)
	
	data_list <- list(n_list = n_list, p_list = p_list, q_list = q_list, k_list = k_list, esR_vector = esR_vector)
	
	return(data_list)
}

# Define a function to apply for each set of parameters 3
plot_evaluate <- function(data_list,rdasim,size,n_iterations){
	
	k_list<-data_list$k_list
	esR_vector <-data_list$esR_vector
	
	h=12
	
	
	countR <- sum(esR_vector == k_list) 
	perR <- 100 * countR / length(k_list) 
	cat("rank is well estimated in :", perR, "% \n")
	
	# Correlation between esR and k
	correlation <- cor(esR_vector, k_list)
	cat("Correlation between esR and k:", signif(correlation, 3), "\n")
	
	perR_vec <- numeric(12)

# チャンクごとに処理
for (i in 1:12) {
  start_idx <- (i - 1) * n_iterations + 1
  end_idx <- i * n_iterations
  
  k_chunk <- k_list[start_idx:end_idx]
  esR_chunk <- esR_vector[start_idx:end_idx]
  
  countR <- sum(esR_chunk == k_chunk)
  perR_e <- 100 * countR / n_iterations
  
  perR_vec[i] <- perR_e
}

	# 
	# #heat
	# 
	# data_table <- table(esR_vector, k_list)
	# full_combination <- expand.grid(esR_vector = 1:h, k_list = 1:10)
	# data_melt <- as.data.frame(as.table(data_table))
	# full_data_melt <- merge(full_combination, data_melt, by = c("esR_vector", "k_list"), all.x = TRUE)
	# full_data_melt$Freq[is.na(full_data_melt$Freq)] <- 0
	# 
	# g<-ggplot(full_data_melt, aes(x = k_list, y = esR_vector, fill = Freq)) +
	# 	geom_tile(color = "black") +
	# 	scale_fill_gradient(low = "white", high = "red", breaks = c(0, 5, 10)) +
	# 	theme_minimal(base_size = 15)  +
	# 	theme(
	# 		axis.title.x = element_text(size = 14),
	# 		axis.title.y = element_text(size = 14),
	# 		plot.title = element_text(hjust = 0, size = 18),
	# 		panel.grid = element_blank(),
	# 		axis.text.x = element_text(size = 10),
	# 		axis.text.y = element_text(size = 10)
	# 	) +
	# 	scale_x_continuous(breaks = 1:10, labels = as.character(1:10)) +
	# 	scale_y_continuous(breaks = 1:h, labels = as.character(1:h))
	# 
	# if (rdasim=="sim1"){
	# 	model = 1
	# } else if (rdasim=="sim2"){
	# 	model = 2
	# } else{
	# 	stop()
	# }
	# 
	# g<-g+labs(title = paste0("Rank Estimate - Model ",model," - Size:",size), x = "True Rank", y = "Estimated Rank", fill = "Frequency") 
	# 
	# print(g) 
	
	return(list(correlation=signif(correlation, 3),perR =perR,perR_vec=perR_vec))
}


tune_sRDA_ust <- function(X, Y, nonzero_vals = c(5, 10, 20, 50, 100), num_vals = c(1, 2, 3), dir_path = ".") {
	set.seed(123)
	tic()
	folds <- createFolds(1:nrow(X), k = 5, list = TRUE)
	results_summary <- data.frame(
		num = integer(),
		nonzero = integer(),
		mean_abs_corr_sum = numeric(),
		mean_MSE = numeric(),
		se_abs_corr_sum = numeric(),
		se_mean_MSE = numeric()
	)
	
	for (num in num_vals) {
		multi <- num > 1
		
		for (nz in nonzero_vals) {
			cat("Running CV for num =", num, ", nonzero =", nz, "\n")
			
			abs_corr_sums <- c()
			mses <- c()
			
			for (i in 1:5) {
				test_idx <- folds[[i]]
				train_idx <- setdiff(1:nrow(X), test_idx)
				
				X_train <- X[train_idx, ]
				Y_train <- Y[train_idx, ]
				X_test <- X[test_idx, ]
				Y_test <- Y[test_idx, ]
				
				# モデル学習
				res <- sRDA(
					predictor = X_train, predicted = Y_train,
					nonzero = nz, penalization = "ust", nr_LVs = num,
					multiple_LV = multi
				)
				
				if (multi) {
					Y_hat <- matrix(0, nrow = nrow(X_test), ncol = ncol(Y_test))
					for (lv in seq_along(res$ALPHA)) {
						xi_test <- scale(X_test %*% res$ALPHA[[lv]])
						Y_hat <- Y_hat + xi_test %*% t(res$BETA[[lv]])
					}
				} else {
					xi_test <- scale(X_test %*% res$ALPHA)
					Y_hat <- xi_test %*% t(res$BETA)
				}
				
				# 評価
				abs_corrs <- sapply(1:ncol(Y_test), function(k) abs(cor(xi_test, Y_test[, k])))
				abs_corr_sum <- sum(abs_corrs)
				mse <- mean((Y_test - Y_hat)^2)
				
				abs_corr_sums <- c(abs_corr_sums, abs_corr_sum)
				mses <- c(mses, mse)
			}
			
			# 結果保存
			results_summary <- rbind(results_summary, data.frame(
				num = num,
				nonzero = nz,
				mean_abs_corr_sum = mean(abs_corr_sums),
				mean_MSE = mean(mses),
				se_abs_corr_sum = sd(abs_corr_sums) / sqrt(length(abs_corr_sums)),
				se_mean_MSE = sd(mses) / sqrt(length(mses))
			))
		}
	}
	
	a <- toc()
	
	# 最良パラメータ取得
	best_row <- results_summary[which.max(results_summary$mean_abs_corr_sum), ]
	best_num <- best_row$num
	best_nonzero <- best_row$nonzero
	
	cat("Best number of LVs:", best_num, "\n")
	cat("Best nonzero value:", best_nonzero, "\n")
	
	saveRDS(results_summary, file.path(dir_path, "ust.rds"))
	saveRDS(a, file.path(dir_path, "time.rds"))
	
	return(list(
		best_num = best_num,
		best_nonzero = best_nonzero,
		time = a
	))
	
}

tune_sRDA_enet <- function(X, Y, 
													 nonzero_vals = c(5, 10, 20, 50, 100), 
													 num_vals = c(1, 2, 3, 4, 5), 
													 lambda_vals = c(0.1, 1, 10, 100),
													 dir_path = ".") {
	set.seed(123)
	tic()
	folds <- createFolds(1:nrow(X), k = 5, list = TRUE)
	
	results_summary <- data.frame(
		num = integer(),
		nonzero = integer(),
		lambda = numeric(),
		mean_abs_corr_sum = numeric(),
		mean_MSE = numeric(),
		se_abs_corr_sum = numeric(),
		se_mean_MSE = numeric()
	)
	
	for (num in num_vals) {
		multi <- num > 1
		
		for (lambda in lambda_vals) {
			for (nz in nonzero_vals) {
				cat("Running CV for num =", num, ", nonzero =", nz, ", lambda =", lambda, "\n")
				
				abs_corr_sums <- c()
				mses <- c()
				
				for (i in 1:5) {
					test_idx <- folds[[i]]
					train_idx <- setdiff(1:nrow(X), test_idx)
					
					X_train <- X[train_idx, ]
					Y_train <- Y[train_idx, ]
					X_test <- X[test_idx, ]
					Y_test <- Y[test_idx, ]
					
					# モデル学習
					res <- sRDA(
						predictor = X_train, predicted = Y_train,
						nonzero = nz, penalization = "enet", nr_LVs = num,
						multiple_LV = multi, ridge_penalty = lambda, cross_validate = FALSE
					)
					
					if (multi) {
						Y_hat <- matrix(0, nrow = nrow(X_test), ncol = ncol(Y_test))
						for (lv in seq_along(res$ALPHA)) {
							xi_test <- scale(X_test %*% res$ALPHA[[lv]])
							Y_hat <- Y_hat + xi_test %*% t(res$BETA[[lv]])
						}
					} else {
						xi_test <- scale(X_test %*% res$ALPHA)
						Y_hat <- xi_test %*% t(res$BETA)
					}
					
					# 評価
					abs_corrs <- sapply(1:ncol(Y_test), function(k) abs(cor(xi_test, Y_test[, k])))
					abs_corr_sum <- sum(abs_corrs)
					mse <- mean((Y_test - Y_hat)^2)
					
					abs_corr_sums <- c(abs_corr_sums, abs_corr_sum)
					mses <- c(mses, mse)
				}
				
				# 結果保存
				results_summary <- rbind(results_summary, data.frame(
					num = num,
					nonzero = nz,
					lambda = lambda,
					mean_abs_corr_sum = mean(abs_corr_sums),
					mean_MSE = mean(mses),
					se_abs_corr_sum = sd(abs_corr_sums) / sqrt(length(abs_corr_sums)),
					se_mean_MSE = sd(mses) / sqrt(length(mses))
				))
			}
		}
	}
	
	a <- toc()
	
	# 最良パラメータ取得
	best_row <- results_summary[which.max(results_summary$mean_abs_corr_sum), ]
	best_num <- best_row$num
	best_nonzero <- best_row$nonzero
	best_lambda <- best_row$lambda
	
	cat("Best number of LVs:", best_num, "\n")
	cat("Best nonzero value:", best_nonzero, "\n")
	cat("Best lambda:", best_lambda, "\n")
	
	saveRDS(results_summary, file.path(dir_path, "enet.rds"))
	saveRDS(a, file.path(dir_path, "enettime.rds"))
	
	return(list(
		best_num = best_num,
		best_nonzero = best_nonzero,
		best_lambda = best_lambda,
		time = a
	))
}



tune_rrda <- function(X, Y, 
											n_lam = 50,
											maxrank = 40,
											dir_path = ".", 
											plot_result = FALSE) {
	set.seed(123)
	tic()
	
	Lambda <- 10^seq(0, 6, length.out = n_lam)
	
	# クロスバリデーション実行
	cv_result <- rrda.cv(Y = Y, X = X, maxrank = maxrank, lambda = Lambda)
	
	# プロット（必要に応じて）
	if (plot_result) {
		rrda.plot(cv_result = cv_result)
	}
	
	a<-toc()
	# 最適パラメータ
	K <- cv_result$opt_min$rank
	L <- cv_result$opt_min$lambda
	mse <- cv_result$opt_min$MSE
	
	cat("Best rank (K):", K, "\n")
	cat("Best lambda (L):", L, "\n")
	cat("Minimum MSE:", mse, "\n")
	
	# 保存
	saveRDS(cv_result, file.path(dir_path, "rrda_cv_result.rds"))
	saveRDS(toc(), file.path(dir_path, "rrda_time.rds"))
	
	# 結果を返す
	return(list(K = K, lambda = L, MSE = mse,time=a))
}

run_parameter_selection <- function(X, Y, split_times = 5, dir_path = ".",nonzero_vals,num_vals,lambda_vals) {
	set.seed(123)  # 再現性のため
	
	rrda_params_list <- list()
	ust_params_list <- list()
	enet_params_list <- list()
	split_folds <- list()
	
	n <- nrow(X)
	
	for (i in 1:split_times) {
		cat("\n=== Split", i, "===\n")
		
		# 10件を検証用にランダムに抽出
		test_idx <- sample(1:n, 10)
		train_idx <- setdiff(1:n, test_idx)
		
		X_train <- X[train_idx, ]
		Y_train <- Y[train_idx, ]
		
		# 保存する fold 情報
		split_folds[[i]] <- list(train_idx = train_idx, test_idx = test_idx)
		
		# ---- RRDA ----
		cat("Tuning RRDA...\n")
		rrda_params <- tune_rrda(X_train, Y_train, dir_path = dir_path, plot_result = FALSE)
		rrda_params_list[[i]] <- rrda_params
		
		# ---- UST ----
		cat("Tuning UST...\n")
		ust_params <- tune_sRDA_ust(X_train, Y_train, dir_path = dir_path,nonzero_vals = nonzero_vals,num_vals = num_vals)
		ust_params_list[[i]] <- ust_params
		
		# ---- ENET ----
		cat("Tuning ENET...\n")
		enet_params <- tune_sRDA_enet(X_train, Y_train, dir_path = dir_path,nonzero_vals=nonzero_vals,num_vals = num_vals,lambda_vals = lambda_vals)
		enet_params_list[[i]] <- enet_params
	}
	
	# 保存
	saveRDS(split_folds, file.path(dir_path, "split_folds.rds"))
	saveRDS(rrda_params_list, file.path(dir_path, "rrda_params_list.rds"))
	saveRDS(ust_params_list, file.path(dir_path, "ust_params_list.rds"))
	saveRDS(enet_params_list, file.path(dir_path, "enet_params_list.rds"))
	
	return(list(
		folds = split_folds,
		rrda = rrda_params_list,
		ust = ust_params_list,
		enet = enet_params_list
	))
}

evaluate_models_mspe <- function(X, Y, results) {
	rrda_mspes <- numeric()
	ust_mspes <- numeric()
	enet_mspes <- numeric()
	
	for (i in seq_along(results$folds)) {
		cat("\n=== Evaluation Split", i, "===\n")
		
		train_idx <- results$folds[[i]]$train_idx
		test_idx <- results$folds[[i]]$test_idx
		
		X_train <- X[train_idx, ]
		Y_train <- Y[train_idx, ]
		X_test <- X[test_idx, ]
		Y_test <- Y[test_idx, ]
		
		### --- RRDA ---
		cat("Evaluating RRDA...\n")
		rrda_param <- results$rrda[[i]]
		res_rrda <- rrda.fit(Y = Y_train, X = X_train,
												 nrank = rrda_param$K, lambda = rrda_param$lambda,
												 scale.X = FALSE, scale.Y = FALSE)
		Y_hat_rrda <- rrda.predict(Bhat = res_rrda, X = X_test)[[1]][[1]][[1]]
		rrda_mspes[i] <- 100 * sum((Y_test - Y_hat_rrda)^2) / (ncol(Y_test) * nrow(Y_test))
		
		### --- UST ---
		cat("Evaluating UST...\n")
		ust_param <- results$ust[[i]]
		multi <- ust_param$best_num > 1
		res_ust <- sRDA(
			predictor = X_train, predicted = Y_train,
			nonzero = ust_param$best_nonzero, penalization = "ust",
			nr_LVs = ust_param$best_num, multiple_LV = multi
		)
		if (multi) {
			Y_hat_ust <- matrix(0, nrow = nrow(X_test), ncol = ncol(Y_test))
			for (lv in seq_along(res_ust$ALPHA)) {
				xi_test <- scale(X_test %*% res_ust$ALPHA[[lv]])
				Y_hat_ust <- Y_hat_ust + xi_test %*% t(res_ust$BETA[[lv]])
			}
		} else {
			xi_test <- scale(X_test %*% res_ust$ALPHA)
			Y_hat_ust <- xi_test %*% t(res_ust$BETA)
		}
		ust_mspes[i] <- 100 * sum((Y_test - Y_hat_ust)^2) / (ncol(Y_test) * nrow(Y_test))
		
		### --- ENET ---
		cat("Evaluating ENET...\n")
		enet_param <- results$enet[[i]]
		multi <- enet_param$best_num > 1
		res_enet <- sRDA(
			predictor = X_train, predicted = Y_train,
			nonzero = enet_param$best_nonzero, penalization = "enet",
			nr_LVs = enet_param$best_num, multiple_LV = multi,
			ridge_penalty = enet_param$best_lambda, cross_validate = FALSE
		)
		if (multi) {
			Y_hat_enet <- matrix(0, nrow = nrow(X_test), ncol = ncol(Y_test))
			for (lv in seq_along(res_enet$ALPHA)) {
				xi_test <- scale(X_test %*% res_enet$ALPHA[[lv]])
				Y_hat_enet <- Y_hat_enet + xi_test %*% t(res_enet$BETA[[lv]])
			}
		} else {
			xi_test <- scale(X_test %*% res_enet$ALPHA)
			Y_hat_enet <- xi_test %*% t(res_enet$BETA)
		}
		enet_mspes[i] <- 100 * sum((Y_test - Y_hat_enet)^2) / (ncol(Y_test) * nrow(Y_test))
	}
	
	# 平均と標準誤差をまとめる
	mspe_summary <- data.frame(
		method = c("rrda", "ust", "enet"),
		mean_mspe = c(mean(rrda_mspes), mean(ust_mspes), mean(enet_mspes)),
		se_mspe = c(
			sd(rrda_mspes) / sqrt(length(rrda_mspes)),
			sd(ust_mspes) / sqrt(length(ust_mspes)),
			sd(enet_mspes) / sqrt(length(enet_mspes))
		)
	)
	
	return(list(
		rrda_mspes = rrda_mspes,
		ust_mspes = ust_mspes,
		enet_mspes = enet_mspes,
		mspe_summary = mspe_summary
	))
}

summarize_times <- function(results) {
	extract_time <- function(time_obj) {
		if (is.null(time_obj)) return(NA)
		return(as.numeric(time_obj$toc - time_obj$tic))
	}
	
	rrda_times <- sapply(results$rrda, function(x) extract_time(x$time))
	ust_times <- sapply(results$ust, function(x) extract_time(x$time))
	enet_times <- sapply(results$enet, function(x) extract_time(x$time))
	
	time_summary <- data.frame(
		method = c("rrda", "ust", "enet"),
		mean_time_sec = c(mean(rrda_times, na.rm = TRUE),
											mean(ust_times, na.rm = TRUE),
											mean(enet_times, na.rm = TRUE)),
		se_time_sec = c(sd(rrda_times, na.rm = TRUE) / sqrt(length(na.omit(rrda_times))),
										sd(ust_times, na.rm = TRUE) / sqrt(length(na.omit(ust_times))),
										sd(enet_times, na.rm = TRUE) / sqrt(length(na.omit(enet_times))))
	)
	
	return(time_summary)
}
