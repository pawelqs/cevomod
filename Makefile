VERSION = `grep -m1 Version DESCRIPTION | cut -f2 -d" "`


release:
	git push
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
	git checkout develop
