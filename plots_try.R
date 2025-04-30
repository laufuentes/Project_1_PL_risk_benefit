library(tidyverse)

colnames(sigma_or) <- paste0("or_",1:ncol(sigma_or))
colnames(sigma_est) <- paste0("est_",1:ncol(sigma_est))

tib <- as_tibble(
    cbind(sigma_or,sigma_est)
) %>%
mutate(id=1:n()) %>% 
pivot_longer(-id, names_to="what", names_pattern="([^_]*)_.*", values_to="values") %>% 
pivot_wider(id_cols=id, names_from="what", values_from="values") %>% 
unnest(cols = c(or, est))

tib %>% 
head(n=1e5)%>%
ggplot()+
geom_point(aes(x=or,y=est),alpha=0.1, color="grey10")+
geom_smooth(aes(x=or,y=est, ))+
geom_abline(intercept=0,slope=1, col="red")

