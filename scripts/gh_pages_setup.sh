#!/bin/bash

echo -e "\033[0;32mSetting up GitHub Pages docs directory\033[0m"

echo -e "Deleting old docs directory"
rm -r docs

echo -e "Moving build/html to docs"
mv build/html docs

echo -e "Adding .nojekyll file to docs (so GitHub pages renders it correctly)"
touch docs/.nojekyll

echo -e "Adding CNAME file"
echo "pyscf.org" >> docs/CNAME

echo -e "Adding google verification to docs"
echo "google-site-verification: google356b16cacc76b077.html" >> docs/google356b16cacc76b077.html
