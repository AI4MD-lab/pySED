"""Small, dependency-free helpers for terminal output."""


def print_section(title, symbol="◆"):
    """Separate stages with one blank line and a shared left margin."""
    gap = " " if symbol == "🚀" else "  "
    print(f"\n    {symbol}{gap}{title}")


def print_detail(label, value):
    """Align parameter names and values across CPU, GPU and plotting output."""
    print(f"        {label:<24}: {value}")


def print_matrix(label, rows):
    """Print cell vectors below their label without changing numeric precision."""
    print(f"        {label}:")
    for row in rows:
        print("            " + "    ".join(f"{value: .8e}" for value in row))


def terminal_label(label):
    """Display Gamma in terminals while retaining math-text labels for plots."""
    return label.replace(r"$\Gamma$", "Γ")
