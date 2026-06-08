import os
import importlib.util

class Translator:
    """
    Handles dynamic localization (l10n).
    Loads requested language strings from locales directory and falls back to en.py if missing.
    """
    def __init__(self, lang_code="en"):
        self.lang_code = str(lang_code).split('_')[0].lower()
        self.strings = self._load_language(self.lang_code)
        
        # Load English as a fallback from locales/en.py if the selected language is not English
        self.fallback_strings = {}
        if self.lang_code != "en":
            self.fallback_strings = self._load_language("en")

    def _load_language(self, code):
        try:
            # Try finding 'locales' relative to this localization.py file
            module_dir = os.path.dirname(os.path.realpath(__file__))
            file_path = os.path.join(module_dir, "locales", f"{code}.py")
            
            # Try finding 'locales' in the current working directory
            if not os.path.exists(file_path):
                cwd_path = os.path.join(os.getcwd(), "locales", f"{code}.py")
                if os.path.exists(cwd_path):
                    file_path = cwd_path
            
            if not os.path.exists(file_path):
                print(f"\nLanguage file not found: {code}.py")
                print(f"Checked paths:\n - {file_path}\n - {cwd_path if 'cwd_path' in locals() else 'N/A'}\n")
                return {}
                
            spec = importlib.util.spec_from_file_location(f"locales_{code}", file_path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module.STRINGS
        except Exception as e:
            print(f"\nFailed to load language file: {e}")
            return {}

    def get(self, key):
        """Returns the localized string. Falls back to English dictionary from en.py or key itself if not found."""
        if self.strings and key in self.strings:
            return self.strings[key]
        elif self.fallback_strings and key in self.fallback_strings:
            return self.fallback_strings[key]
        else:
            return f"[{key}]"
