using System.Runtime.InteropServices;
using System.Windows;

namespace EllipticCurves.Explorer.Controls;

internal static class ClipboardActions
{
    public static bool CopyText(Window? owner, string text, string title)
    {
        try
        {
            Clipboard.SetText(text);
            return true;
        }
        catch (ExternalException)
        {
            ConfirmationWindow.ShowMessage(owner, title, "The clipboard is busy. Please try again.");
            return false;
        }
    }
}
