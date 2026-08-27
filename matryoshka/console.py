"""Rich terminal presentation for the public analysis workflow."""

from __future__ import annotations

from collections.abc import Callable
from contextlib import AbstractContextManager
from pathlib import Path
from types import TracebackType

from rich.console import Console
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

from .detect import DetectorRuntime
from .detector_workflow import DetectorExecution


class ReferenceProgress(AbstractContextManager[Callable[[str, int, int], None]]):
    """Adapt the reference scanner callback to a Rich progress bar."""

    def __init__(self, console: Console, profile: str, threads: int, quiet: bool) -> None:
        self._progress = Progress(
            SpinnerColumn(),
            TextColumn("[bold cyan]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TimeElapsedColumn(),
            console=console,
            disable=quiet,
        )
        self._task = self._progress.add_task(
            f"Scanning {profile} component references ({threads} workers)",
            total=None,
        )

    def __enter__(self) -> Callable[[str, int, int], None]:
        self._progress.start()
        return self.update

    def update(self, reference: str, completed: int, total: int) -> None:
        self._progress.update(
            self._task,
            completed=completed,
            total=total,
            description=f"Scanning references · {reference}",
        )

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self._progress.stop()


class WorkflowConsole:
    """Small presentation layer that keeps scientific work out of the CLI."""

    def __init__(self, *, quiet: bool = False) -> None:
        self.quiet = quiet
        self.console = Console(stderr=True, highlight=False)

    def header(
        self,
        fasta: Path,
        out: Path,
        profile: str,
        threads: int,
        detector_mode: str,
    ) -> None:
        if self.quiet:
            return
        self.console.print(Panel.fit(
            "\n".join([
                f"[bold]Input[/]      {fasta}",
                f"[bold]Results[/]    {out}",
                f"[bold]Rules[/]      {profile}",
                f"[bold]Detectors[/]  {detector_mode}",
                f"[bold]Threads[/]    {threads}",
            ]),
            title="[bold magenta]Matryoshka analysis[/]",
            border_style="magenta",
        ))

    def stage(self, number: int, total: int, title: str) -> None:
        if not self.quiet:
            self.console.rule(f"[bold magenta]{number}/{total}[/] [bold]{title}[/]")

    def info(self, message: str) -> None:
        if not self.quiet:
            self.console.print(f"[cyan]•[/] {message}")

    def warning(self, message: str) -> None:
        if not self.quiet:
            self.console.print(f"[yellow]![/] {message}")

    def success(self, message: str) -> None:
        if not self.quiet:
            self.console.print(f"[green]✓[/] {message}")

    def reference_progress(self, profile: str, threads: int) -> ReferenceProgress:
        return ReferenceProgress(self.console, profile, threads, self.quiet)

    def detector_summary(self, records: list[DetectorExecution]) -> None:
        if self.quiet:
            return
        table = Table(title="Optional detector status", box=None, pad_edge=False)
        table.add_column("Detector", style="bold")
        table.add_column("Status")
        table.add_column("Version")
        table.add_column("Output / reason", overflow="fold")
        colours = {
            "completed": "green",
            "provided": "cyan",
            "disabled": "dim",
            "unavailable": "yellow",
            "failed": "red",
        }
        for record in records:
            colour = colours.get(record.status, "white")
            version = record.tool_version or "—"
            if record.database_version:
                version += f" · DB {record.database_version}"
            detail = record.output or record.message or record.runtime
            table.add_row(
                record.label,
                f"[{colour}]{record.status}[/{colour}]",
                version,
                detail,
            )
        self.console.print(table)

    def preflight_summary(
        self,
        core_ready: bool,
        blastn_path: str | None,
        runtimes: list[tuple[str, str, DetectorRuntime]],
    ) -> None:
        """Show the core and optional-tool readiness in one readable table."""
        table = Table(title="Matryoshka preflight", pad_edge=False)
        table.add_column("Component")
        table.add_column("Status")
        table.add_column("Runtime / reason", overflow="fold")
        core_status = "[green]ready[/]" if core_ready else "[red]missing[/]"
        table.add_row("BLAST+ component scan", core_status, blastn_path or "blastn not found")
        for _, label, runtime in runtimes:
            status = "[green]available[/]" if runtime.available else "[yellow]optional[/]"
            detail = (
                runtime.executable
                if runtime.source == "path"
                else f"{runtime.source}: {runtime.reason}"
            )
            table.add_row(label, status, detail)
        self.console.print(table)
        self.console.print(
            "[dim]Detector modes: available = best effort · all = require every tool · "
            "none = component/reference rules only[/]"
        )

    def contig(self, name: str, length: int, roots: int, loci: int) -> None:
        self.success(
            f"{name}: {length:,} bp · {roots} top-level annotations · {loci} locus views"
        )

    def finish(self, features: int, loci: int, records: int, out: Path) -> None:
        if self.quiet:
            return
        table = Table.grid(padding=(0, 2))
        table.add_column(style="bold")
        table.add_column(justify="right")
        table.add_row("FASTA records", str(records))
        table.add_row("Annotated features", str(features))
        table.add_row("Locus maps and tables", str(loci))
        table.add_row("Results directory", str(out))
        self.console.print(Panel(
            table,
            title="[bold green]Analysis complete[/]",
            border_style="green",
        ))
