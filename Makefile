VERSION = `grep -m1 Version DESCRIPTION | cut -f2 -d" "`


release:
	git push
	R CMD INSTALL --preclean --no-multiarch --with-keep.source ../cevomod
	git checkout master
	git pull
	git merge develop
	git tag -a $(VERSION) -m "version $(VERSION)"
	git push --tags
	git checkout develop


merge-to-main:
	git push
	git checkout master
	git pull
	git merge develop
	git push
	git checkout develop


run-cevomod-on-cevodatasets:
	Rscript inst/run_cevomod_on_cevodatasets.R


run-browser:
	R -e 'cevomod::run_browser()'

