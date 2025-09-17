param (
    [string]$log_path
)

[float]$sample_delay = 1
[float]$unit_factor = 1/(1024*1024)
$log_file = Join-Path -Path $log_path -ChildPath 'memory.txt'
$stop_file =Join-Path -Path $log_path -ChildPath 'memory.stop'

$title = 'Memory Profiler'
$Host.UI.RawUI.WindowTitle = $title


Write-Output $pid
while ($true)
{
    $memory = (Get-CimInstance -ClassName Win32_OperatingSystem).FreePhysicalMemory*$unit_factor
    $memory | Add-Content -Path $log_file
    '{0:N2} GB free' -f $memory | Write-Output
    if (Test-Path -Path $stop_file) {break}
    Start-Sleep -Seconds $sample_delay
}