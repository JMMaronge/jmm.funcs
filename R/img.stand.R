#' A function to standardize NIfTI objects
#'
#' This function allows you to standardize NIfTI objects. It also allows for you to do a trimmed standardizating (subtract trimmed mean, divide by trimmed standard deviation)
#' @param image NIfTI object, image to standardize.
#' @param mask NIfTI object, brain mask of image.
#' @param trim scalar between (0,1), percent of trimming to do from the highest and lowest intensities
#' @keywords NIfTI, normalization, MRI, trimmed
#' @export
#' @examples
#' img.stand()


img.stand<-function(image,mask,trim){
  require(fslr)
  require(DescTools)
  dat<-data.frame(flair=c(image),b.mask=c(mask))
  trimmed.dat<-Trim(dat$flair[dat$b.mask==1],trim)
  trimmed.mean<-mean(trimmed.dat)
  trimmed.sd<-sd(trimmed.dat)
  norm.img<-(image-trimmed.mean)/trimmed.sd
  return(norm.img)
}