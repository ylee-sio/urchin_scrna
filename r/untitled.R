t1 = read_csv("~/singapore_data.csv")
t2 = mutate(t1, up_down = case_when(log2foldchange > 0 ~ "UP",log2foldchange < 0 ~ "DOWN", -log10(padjust) < 1.3 ~ "NO"))

p = ggplot(t2) +
geom_point(aes(x=log2foldchange, y=padjust, col = up_down)) +
geom_hline(yintercept = 1.3) +
geom_vline(xintercept = 0)