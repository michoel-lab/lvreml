from setuptools import setup,find_packages

setup(name='lvreml',
	version='0.1',
	description='Python implementation of LVREML Package',
	url='https://github.com/michoel-lab/lvreml',
	author='Tom Michoel',
	author_email='tom.michoel@uib.no',
	license='MIT',
	packages=['lvreml', 'lvreml/modules','lvreml/figcodes'],
	zip_safe=False)
