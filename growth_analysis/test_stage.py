"""Synthetic integrity and dispatcher checks; no saved scientific output replay."""
from pathlib import Path
import hashlib,json,tempfile,unittest
from unittest.mock import patch
from growth_analysis.manifest import verify_manifest


class StageTests(unittest.TestCase):
    def test_manifest_rejects_changed_source(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);data=root/'tiny.bin';data.write_bytes(b'frozen synthetic bytes')
            manifest=root/'manifest.json';manifest.write_text(json.dumps({'files':[dict(path='tiny.bin',absolute=False,bytes=data.stat().st_size,sha256=hashlib.sha256(data.read_bytes()).hexdigest())]}))
            self.assertEqual(verify_manifest(manifest,root)['verified_files'],1)
            data.write_bytes(b'changed synthetic bytes')
            with self.assertRaises(ValueError):verify_manifest(manifest,root)

    def test_fixed_entrypoint_dispatches_analysis_without_training(self):
        import study.run as fixed
        with patch('growth_analysis.analyze.main') as analyze:
            fixed.main()
            analyze.assert_called_once_with(Path(fixed.__file__).resolve().parents[1]/'config.json')


if __name__=='__main__':unittest.main()
