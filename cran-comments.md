## Test environments
* Ubuntu 16.04 (travis): r-devel (2019-09-18 r77193), r-release, r-oldrelease
* winbuilder: r-release (3.6.1)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:
Maintainer: ‘Alexander Rix <alexrix@umich.edu>’

This is actually an update and *name change request* to a package I currently
maintain on CRAN, https://cran.r-project.org/package=hdbm.

My understanding from browsing the R mailing list is that package name changes
are only granted for a good reason. The reason I am requesting the name change
is that my PI doesn't like the name 'hdbm'. I asked her to confirm the name
weeks ago when the 'hdbm' was initially submitted, but she apparently doesn't
always read my emails.

I doubt this is a good reason from CRAN's perspective, but I think my only
alternative is to submit 'bama' as a new package, and then to deprecate 'hdbm'
and point it to 'bama'. This is not entirely satisfactory to me, as the name
'hdbm' will be unavailable for others to use.
