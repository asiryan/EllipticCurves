using System.ComponentModel;
using System.Diagnostics;
using System.Windows;

namespace EllipticCurves.Explorer.Controls;

internal static class BrowserActions
{
    public static void Open(Window owner, string url, string title)
    {
        try
        {
            using var browser = Process.Start(new ProcessStartInfo(url) { UseShellExecute = true });
        }
        catch (Exception error) when (error is Win32Exception or InvalidOperationException)
        {
            ConfirmationWindow.ShowMessage(owner, title,
                "Could not open the browser. Open this address manually:\n" + url);
        }
    }
}
