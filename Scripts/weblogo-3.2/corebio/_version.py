

# Keywords between dollar signs are substituted by subversion.
# The date and build will only tell the truth after a branch or tag,
# since different files in trunk will have been changed at different times
date ="$Date: 2012-01-30 17:55:21 -0800 (Mon, 30 Jan 2012) $".split()[1]
revision = "$Revision: 129 $".split()[1]


__version__ = '3.2' 


description = "CoreBio %s (%s)" % (__version__,  date)


