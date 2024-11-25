$VER = "v0.02.000"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* Complete rewrite
* NEW: interface takes a results workbook and asks for the Stage3.csv file from the experiment
* NEW: generates a 384-well custom manifest that can be loaded directly into the main pipeline for reruns
"@

Move-Item -Path "changelog.txt" -Destination "changelog_old.txt"
Set-Content -Path "changelog.txt" -Value $VER
Add-Content -Path "changelog.txt" -Value $DATE
Add-Content -Path "changelog.txt" -Value $COMMENT
Add-Content -Path "changelog.txt" -Value ""
Get-Content -Path "changelog_old.txt" | Add-Content -Path "changelog.txt"
Remove-Item -Path "changelog_old.txt"

git add -u
$MSG = $COMMENT
git commit -m $MSG
git tag -a $VER -m $MSG

git push
git push origin $VER
