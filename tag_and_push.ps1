$VER = "v0.03.002"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* FIXED: column counting which stopped successful parsing of workbooks
* NEW: added a clear screen button to refresh the state of the application
* CHANGED: rebuilt the uploaders to use callbacks for better clarity
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
