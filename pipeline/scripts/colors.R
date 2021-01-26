library(colorspace)

divergingValueToColor <- function( x, lim=2, gran=10  )
{
  origo = gran+1
  # generate pre-tested Blue-Gray-Red palette
  colcount = gran*2+1
  pal <- diverging_hcl(n=colcount,h1=255,h2=12,c1=130,
                       cmax=80,l1=31,l2=78,p1=1,p2=1.3)
  
  # set extreme thresholds (values must fall between these two)
  x[ x > lim ] = lim
  x[ x < (-lim) ] = -lim
  
  # calculate stepping
  step <- lim/gran
  
  # lookup colors from paletteand return this vector
  return( pal[ origo+round(x/step) ] )
}

isotype_cols <- c( IgD=hcl( 359,1,50 ),
                   IgM=hcl( 60,0,80 ),
                   IgG=hcl( 30,97,99 ),
                   IgA=hcl( 41,60,99 ),
                   IgE=hcl( 60,98,50 ) )