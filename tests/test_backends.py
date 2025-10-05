# tests/test_backends.py
import unittest
from unittest.mock import patch, Mock, MagicMock
import subprocess
from ..backends.schedulers.local import LocalScheduler
from ..backends.schedulers.slurm import SlurmScheduler
from ..backends.launchers.srun import SrunLauncher
from ..backends.launchers.mpirun import MpiRunLauncher
from ..core.types import JobID, JobStatus

class TestLocalScheduler(unittest.TestCase):
    def setUp(self):
        self.scheduler = LocalScheduler()

    def test_submit(self):
        job_script = MagicMock(spec=Path)
        job_id = self.scheduler.submit(job_script)
        self.assertTrue(str(job_id).startswith('local_'))
        self.assertIsInstance(job_id, JobID)

    @patch('subprocess.run')
    def test_wait_success(self, mock_run):
        job_script = MagicMock(spec=Path)
        mock_run.return_value = Mock(returncode=0, stdout='success', stderr='ok')
        job_id = self.scheduler.submit(job_script)
        result = self.scheduler.wait(job_id)
        self.assertEqual(result['status'], JobStatus.COMPLETED)
        self.assertEqual(result['stdout'], 'success')

    @patch('subprocess.run')
    def test_wait_timeout(self, mock_run):
        job_script = MagicMock(spec=Path)
        mock_run.side_effect = subprocess.TimeoutExpired
        job_id = self.scheduler.submit(job_script)
        result = self.scheduler.wait(job_id, timeout=1)
        self.assertEqual(result['status'], JobStatus.FAILED)
        self.assertIn('Timeout', result['error'])

    def test_monitor(self):
        job_id = JobID('test')
        status = self.scheduler.monitor(job_id)
        self.assertEqual(status, JobStatus.COMPLETED)

    def test_cancel(self):
        job_id = JobID('test')
        success = self.scheduler.cancel(job_id)
        self.assertTrue(success)

class TestSlurmScheduler(unittest.TestCase):
    @patch('subprocess.run')
    def test_submit_success(self, mock_run):
        mock_run.return_value = Mock(returncode=0, stdout='Submitted batch job 123')
        job_script = MagicMock(spec=Path)
        scheduler = SlurmScheduler()
        job_id = scheduler.submit(job_script)
        self.assertEqual(job_id, JobID('123'))
        mock_run.assert_called_once_with(['sbatch', str(job_script)], capture_output=True, text=True, check=True, env=None)

    @patch('subprocess.run')
    def test_monitor_pending(self, mock_run):
        mock_run.side_effect = [
            Mock(returncode=0, stdout='PD\n'),  # squeue
            Mock(returncode=0, stdout='')  # sacct not called
        ]
        scheduler = SlurmScheduler()
        status = scheduler.monitor(JobID('123'))
        self.assertEqual(status, JobStatus.PENDING)
        self.assertEqual(mock_run.call_count, 1)

    @patch('subprocess.run')
    def test_monitor_completed(self, mock_run):
        mock_run.side_effect = [
            Mock(returncode=1, stdout=''),  # squeue no job
            Mock(returncode=0, stdout='COMPLETED')  # sacct
        ]
        scheduler = SlurmScheduler()
        status = scheduler.monitor(JobID('123'))
        self.assertEqual(status, JobStatus.COMPLETED)

    @patch('subprocess.run')
    def test_cancel_success(self, mock_run):
        mock_run.return_value = Mock(returncode=0)
        scheduler = SlurmScheduler()
        success = scheduler.cancel(JobID('123'))
        self.assertTrue(success)
        mock_run.assert_called_once_with(['scancel', '123'], capture_output=True, check=True, env=None)

    @patch('subprocess.run')
    def test_wait_completed(self, mock_run):
        # Simulate polling: first pending, then completed
        mock_run.side_effect = [
            Mock(returncode=0, stdout='PD'),  # monitor pending
            Mock(returncode=1, stdout=''),    # squeue empty
            Mock(returncode=0, stdout='COMPLETED'),  # sacct
            Mock(returncode=0, stdout='JobID|...|COMPLETED|...')  # sacct details
        ]
        scheduler = SlurmScheduler()
        job_id = JobID('123')
        # Submit mock
        with patch.object(scheduler, 'submit', return_value=job_id):
            result = scheduler.wait(job_id)
        self.assertEqual(result['status'], JobStatus.COMPLETED)
        self.assertIn('sacct_details', result)

class TestSrunLauncher(unittest.TestCase):
    @patch('subprocess.run')
    def test_launch_success(self, mock_run):
        mock_run.return_value = Mock(returncode=0, stdout='output', stderr='3.45')
        launcher = SrunLauncher()
        result = launcher.launch(['app'], nodes=1, procs_per_node=4)
        self.assertEqual(result['runtime'], 3.45)
        self.assertEqual(result['returncode'], 0)

    @patch('subprocess.run')
    def test_launch_failed(self, mock_run):
        mock_run.return_value = Mock(returncode=1, stdout='err out', stderr='err')
        launcher = SrunLauncher()
        result = launcher.launch(['app'], nodes=1, procs_per_node=4)
        self.assertEqual(result['runtime'], 0.0)
        self.assertEqual(result['status'], 'FAILED')

class TestMpiRunLauncher(unittest.TestCase):
    @patch('subprocess.run')
    def test_launch_success(self, mock_run):
        mock_run.return_value = Mock(returncode=0, stdout='output', stderr='real 5.67')
        launcher = MpiRunLauncher()
        result = launcher.launch(['app'], nodes=1, procs_per_node=4)
        self.assertEqual(result['runtime'], 5.67)
        self.assertEqual(result['returncode'], 0)

    @patch('subprocess.run')
    def test_launch_failed(self, mock_run):
        mock_run.return_value = Mock(returncode=1, stdout='err out', stderr='err')
        launcher = MpiRunLauncher()
        result = launcher.launch(['app'], nodes=1, procs_per_node=4)
        self.assertEqual(result['runtime'], 0.0)
        self.assertEqual(result['status'], 'FAILED')

if __name__ == '__main__':
    unittest.main()
