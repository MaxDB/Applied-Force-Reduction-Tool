param (
    [string]$log_path,
    [int16]$display = 1
)


[float]$sample_delay = 1
[float]$unit_factor = 1/(1024*1024)
$log_file = Join-Path -Path $log_path -ChildPath 'memory.txt'
$stop_file =Join-Path -Path $log_path -ChildPath 'memory.stop'

$title = 'Memory Profiler'
$Host.UI.RawUI.WindowTitle = $title

$start_time = Get-Date
if ($display -eq 1) {$start_time | Write-Output}
$start_time | Add-Content -Path $log_file
$sample_delay | Add-Content -Path $log_file
"--" | Add-Content -Path $log_file
while ($true)
{
    $memory = (Get-CimInstance -ClassName Win32_OperatingSystem).FreePhysicalMemory*$unit_factor
    $memory | Add-Content -Path $log_file
    if ($display -eq 1) {'{0:N2} GB free' -f $memory | Write-Output}
    if (Test-Path -Path $stop_file) {break}
    Start-Sleep -Seconds $sample_delay
}

$stop_time = Get-Date
"--" | Add-Content -Path $log_file
$stop_time | Add-Content -Path $log_file

if ($display -eq 1){
    $stop_time | Write-Output
    Write-Host "Press any key to exit..."
    $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown") | Out-Null
}