#' Split Rows
#' 
#' Split rows in data.frame if split point falls between start and end of range in data.frame
#' @param df A data.frame
#' @param splitAt A number specifing the point within range to split at
#' @param id The name of id column, default is "id"
#' @param start The name of start of range column, default is "start"
#' @param end The name of end of range column, default is "end"
#' @return A data.frame with new start and end range, group, start_type and end_type
#' @author Johan Junkka
#' @export

split_rows <- function(df, splitAt, id ="id", start = "start", end = "end"){
  
  if (!is.data.frame(df)) stop("df must be a data.frame")

  # Determin if row should be split
  splitRow <- ifelse( df[ ,start] < splitAt & df[ ,end] > splitAt, TRUE, FALSE)
  splitRow[is.na(splitRow)] <- FALSE
  
  # for each row check if split and make new start
  dfid      <- df[ ,id]
  df.end  <- df[ ,end]
  
  rowOne <- data.frame(
    id         = dfid,
    start      = df[ ,start],
    end        = ifelse(splitRow == TRUE, splitAt, df.end),
    group      = ifelse(df[ ,start] < splitAt, 0, splitAt),
    start_type = rep(0, length(df[ ,id])),
    end_type   = ifelse(splitRow == TRUE, 1, 0)
  )
  
  # make new
  new.id <- c()
  new.end <- c()
  
  for (j in 1:length(dfid)){
    # if splitRow is true populate new vectors
    if (splitRow[j] == TRUE) {
      k <- length(new.id)+1
      new.id[k]  <- dfid[j]
      new.end[k] <- df.end[j]
    }
  }
  
  rowTwo <- data.frame(
    id          = new.id,
    start       = rep(splitAt, length(new.id)),
    end         = new.end,
    group       = rep(splitAt, length(new.id)),
    start_type  = rep(1, length(new.id)),
    end_type    = rep(0, length(new.id))
  )
  
  df <- rbind(rowOne,rowTwo)
  
  return(df[order(df$id, df$start), ])
  
}


#' Age split event history data
#'
#' Split event history data into a sequence of age groups. A row is split if it the range 
#'   between start and end stretches over a age group. Each row is assigned a new start 
#'   and end point, a group and a duration within that group.
#' @param dat data.frame containing event history data
#' @param  start = "start" Name of start column
#' @param  end = "end" Name of end column
#' @param  id = "id" Name of id column
#' @param  event = "event" Name of event column
#' @param  split_by = 5 Age group size
#' @export

age_split <- function(dat, start = "start", end = "end", id = "id", event = "event", split_by = 5) {

  id    <- dat[ ,id]
  start <- dat[ ,start]
  end   <- dat[ ,end]
  event <- dat[ ,event]

  # for each range calculate number of split points

  # first split
  first     <- ((start%/%split_by)+1)*split_by
  # if split at all
  splitRow  <- ifelse( end <= first, FALSE, TRUE)


  start2    <- start[splitRow]
  end2      <- end[splitRow]
  id2       <- id[splitRow]
  first2    <- first[splitRow]
  event2    <- event[splitRow]

  no_split <- data.frame(
    id    = id[!splitRow],
    start = start[!splitRow],
    end   = end[!splitRow],
    group = first[!splitRow] - split_by,
    event = event[!splitRow]
  )

  if (all(!splitRow)) {
    warning('No range falls within split_by.', call. = FALSE) 
    dat$id            <- id
    no_split$duration <- no_split$end - no_split$start
    no_split          <- join(no_split, dat, by="id", type="left")
    no_split          <- no_split[no_split$duration > 0, ]

    return(no_split[order(no_split$id, no_split$start), ])
  }
  new_start <- c()  
  new_end   <- c()
  new_id    <- c()
  group     <- c()
  new_event <- c()
  k         <- 0

  for (j in 1:length(start2)) {
    # new vector key
    k <- k + 1
    # push first 
    new_start[k]  <- start2[j]
    new_end[k]    <- first2[j]
    new_id[k]     <- id2[j]
    group[k]      <- first2[j] - split_by # split_by iteration
    new_event[k]  <- 0
    # how many splits
    range  <- end2[j] - first2[j]
    splits <- range %/% split_by

    # iterate over splits and push into new vector
    if (splits > 0) {
      for (i in 1:splits){
        k <- k + 1
        # new start is first2 split + split_by * split oreder -1
        new_start[k]  <- first2[j] + (split_by * (i - 1))
        new_end[k]    <- first2[j] + (split_by * i)
        new_id[k]     <- id2[j]
        group[k]      <- new_start[k] 
        new_event[k]  <- 0
      }
    }
    k <- k + 1
    # push last
    new_start[k]  <- first2[j] + (split_by * splits)
    new_end[k]    <- end2[j]
    new_id[k]     <- id2[j]
    group[k]      <- new_start[k] 
    new_event[k]  <- event2[j]
  }

  split <- data.frame(
    id    = new_id,
    start = new_start,
    end   = new_end,
    group = group,
    event = new_event
  )

  out <- rbind(no_split,split)
  out$duration <- out$end - out$start
  dat$id <- id
  out <- plyr::join(out, dat, by="id", type="left")
  out <- out[out$duration > 0, ]

  return(out[order(out$id, out$start), ])

}