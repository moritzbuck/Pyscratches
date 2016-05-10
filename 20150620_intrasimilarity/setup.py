from distutils.core import setup

setup(
      name             = 'intrasimilarity',
      version          = '0.0',
      description      = 'assembly comparison tools',
      long_description = 'an other day',
      license          = 'MIT',
      url              = 'https://github.com/viklund/AssemblyValidation',
      author           = 'Moritz Buck',
      author_email     = 'mrtz.buck@gmail.com',
      classifiers      = ['Topic :: Scientific/Engineering :: Bioinformatics'],
      packages         = ['intrasimilarity'],
      requires         = ['matplotlib', 'scipy', 'pandas', 'Bio'],
)
