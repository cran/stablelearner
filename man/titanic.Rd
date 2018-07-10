\name{titanic}
\alias{titanic}

\title{Passengers and Crew on the RMS Titanic}

\description{
  the \code{titanic} data is a complete list of passengers and crew members on 
  the RMS Titanic. It includes a variable indicating whether a person did 
  survive the sinking of the RMS Titanic on April 15, 1912.
}

\usage{data("titanic")}

\format{
A data frame containing 2207 observations on 11 variables.
\describe{
  \item{name}{a string with the name of the passenger.} 
  \item{gender}{a factor with levels \code{male} and \code{female}.} 
  \item{age}{a numeric value with the persons age on the day of the 
  sinking. The age of babies (under 12 months) is given as a fraction of one 
  year (1/month).} 
  \item{class}{a factor specifying the class for passengers or the type of 
  service aboard for crew members.}
  \item{embarked}{a factor with the persons place of of embarkment.}
  \item{country}{a factor with the persons home country.}
  \item{ticketno}{a numeric value specifying the persons ticket number 
  (\code{NA} for crew members).}
  \item{fare}{a numeric value with the ticket price (\code{NA} for crew 
  members, musicians and employees of the shipyard company).}
  \item{sibsp}{an ordered factor specifying the number if siblings/spouses 
  aboard; adopted from Vanderbild data set (see below).}
  \item{parch}{an ordered factor specifying the number of parents/children 
  aboard; adopted from Vanderbild data set (see below).}
  \item{survived}{a factor with two levels (\code{no} and \code{yes}) 
  specifying whether the person has survived the sinking.}
}
}

\details{

The website \url{http://www.encyclopedia-titanica.org/} offers detailed 
information about passengers and crew members on the RMS Titanic. According to
the website 1317 passengers and 890 crew member were abord.

8 musicians and 9 employees of the shipyard company are listed as 
passengers, but travelled with a free ticket, which is why they have \code{NA}
values in \code{fare}. In addition to that, \code{fare} is truely missing for 
a few regular passengers.

}

\source{
  The complete list of persons on the RMS titanic was downloaded from 
  \url{http://www.encyclopedia-titanica.org/} on April 5, 2016. The 
  information given in \code{sibsp} and \code{parch} was adopoted from a data
  set obtained from \url{http://biostat.mc.vanderbilt.edu/DataSets}.
}

\references{
  \url{http://www.encyclopedia-titanica.org/} and 
  \url{http://biostat.mc.vanderbilt.edu/DataSets}.
}

\examples{
data("titanic", package = "stablelearner")
summary(titanic)
}

\keyword{datasets}
