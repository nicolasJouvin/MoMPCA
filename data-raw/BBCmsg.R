## code to prepare `BBCmsg` dataset goes here
preprocess_corpus = function(corpus, stop_words='english', language='english'){
  corpus <- tm::tm_map(corpus, tm::content_transformer(tolower))
  corpus <- tm::tm_map(corpus, tm::content_transformer(removeNumbers))
  #corpus <- tm_map(corpus, tm::content_transformer(removeWords),
  #                 stop_words)

  corpus <- tm::tm_map(corpus, tm::removePunctuation)
  corpus <- tm::tm_map(corpus, tm::stemDocument, language = language)
  corpus <- tm::tm_map(corpus, tm::stripWhitespace)
  return(corpus)
}

msg1 = as.character(read.table('msgA.txt',stringsAsFactors = FALSE))
msg2 = as.character(read.table('msgB.txt',stringsAsFactors = FALSE))
msg3 = as.character(read.table('msgC.txt',stringsAsFactors = FALSE))
msg4 = as.character(read.table('msgD.txt',stringsAsFactors = FALSE))
corpus = list()
for (k in 1:4) {
  do.call("<-",list('msg', as.name(paste0('msg', k))))
  corpus[[k]] = paste(msg, collapse =  ' ')
}

corpus = tm::Corpus(tm::VectorSource(corpus))
corpus = preprocess_corpus(corpus = corpus, language = 'english')
corpus = tm::tm_map(corpus, tm::content_transformer(removeWords), tm::stopwords('english'))
corpus = tm:tm_map(corpus, tm::stripWhitespace)

allmsg = c()
for (k in 1:4) {
  allmsg = c( allmsg, unlist(strsplit(corpus[[k]]$content, ' ')))
  msg = unlist(strsplit(corpus[[k]]$content, ' '))
  do.call("<-",list( as.name(paste0('msg', k)), as.name('msg')))
}
n1 = n2 = 0;
for (k in 1:4) {
  do.call("<-",list('msg', as.name(paste0('msg', k))))
  n1 = n1 + length(msg)
}
stopifnot(n1 == length(allmsg))

BBCmsg = list(msg1 = msg1, msg2 = msg2, msg3 = msg3, msg4 = msg4)


usethis::use_data("BBCmsg")
