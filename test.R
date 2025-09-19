#test
library(wordcloud2)
library(webshot)
webshot::install_phantomjs()
library("htmlwidgets")

my_graph <- wordcloud2(data=read_delim("test_wordcloud.csv",";",show_col_types = F),
                       size=1.5, color = c("#4A7EBB","#7F7F7F","#4A7EBB",
                                           "#72A16E",
                                           "#FF0000", "#FF0000", "#FF0000"),
                       minRotation = 0, maxRotation = 0,rotateRatio = 0)


my_graph <- wordcloud2(data=read_delim("test_wordcloud.csv",";",show_col_types = F),
           size = 1, minRotation = -pi/6, maxRotation = -pi/6, rotateRatio = 1)

saveWidget(my_graph,"/Users/piefouca/Desktop/tmp.html",selfcontained = F)
webshot("/Users/piefouca/Desktop/tmp.html","/Users/piefouca/Desktop/test_wordcloud.png",
        delay =5, vwidth = 480, vheight=480)
