class Solution {
 public:
  vector<int> twoSum(vector<int>& nums, int target) {
    unordered_map<int, int> numToIndex;

    for (int i = 0; i < nums.size(); ++i) {
      if (const auto it = numToIndex.find(target - nums[i]);
          it != numToIndex.cend())
        return {it->second, i};
      numToIndex[nums[i]] = i;
    }

    throw;
  }
};
class Solution {
 public:
  bool isPalindrome(int x) {
    if (x < 0)
      return false;

    long reversed = 0;
    int y = x;

    while (y > 0) {
      reversed = reversed * 10 + y % 10;
      y /= 10;
    }

    return reversed == x;
  }
};
class Solution {
 public:
  int romanToInt(string s) {
    int ans = 0;
    vector<int> roman(128);

    roman['I'] = 1;c
    roman['V'] = 5;
    roman['X'] = 10;
    roman['L'] = 50;
    roman['C'] = 100;
    roman['D'] = 500;
    roman['M'] = 1000;

    for (int i = 0; i + 1 < s.length(); ++i)
      if (roman[s[i]] < roman[s[i + 1]])
        ans -= roman[s[i]];
      else
        ans += roman[s[i]];

    return ans + roman[s.back()];
  }
};
class Solution {
 public:
  string longestCommonPrefix(vector<string>& strs) {
    if (strs.empty())
      return "";

    for (int i = 0; i < strs[0].length(); ++i)
      for (int j = 1; j < strs.size(); ++j)
        if (i == strs[j].length() || strs[j][i] != strs[0][i])
          return strs[0].substr(0, i);

    return strs[0];
  }
};
class Solution {
 public:
  bool isValid(string s) {
    stack<char> stack;

    for (const char c : s)
      if (c == '(')
        stack.push(')');
      else if (c == '{')
        stack.push('}');
      else if (c == '[')
        stack.push(']');
      else if (stack.empty() || pop(stack) != c)
        return false;

    return stack.empty();
  }

 private:
  int pop(stack<char>& stack) {
    const int c = stack.top();
    stack.pop();
    return c;
  }
};
class Solution {
 public:
  ListNode* mergeTwoLists(ListNode* list1, ListNode* list2) {
    if (!list1 || !list2)
      return list1 ? list1 : list2;
    if (list1->val > list2->val)
      swap(list1, list2);
    list1->next = mergeTwoLists(list1->next, list2);
    return list1;
  }
};
class Solution {
 public:
  ListNode* mergeTwoLists(ListNode* list1, ListNode* list2) {
    if (!list1 || !list2)
      return list1 ? list1 : list2;
    if (list1->val > list2->val)
      swap(list1, list2);
    list1->next = mergeTwoLists(list1->next, list2);
    return list1;
  }
};
class Solution {
 public:
  int removeDuplicates(vector<int>& nums) {
    int i = 0;

    for (const int num : nums)
      if (i < 1 || num > nums[i - 1])
        nums[i++] = num;

    return i;
  }
};
class Solution {
 public:
  int removeElement(vector<int>& nums, int val) {
    int i = 0;

    for (const int num : nums)
      if (num != val)
        nums[i++] = num;

    return i;
  }
};
class Solution {
 public:
  int strStr(string haystack, string needle) {
    const int m = haystack.length();
    const int n = needle.length();

    for (int i = 0; i < m - n + 1; i++)
      if (haystack.substr(i, n) == needle)
        return i;

    return -1;
  }
};
class Solution {
 public:
  int searchInsert(vector<int>& nums, int target) {
    int l = 0;
    int r = nums.size();

    while (l < r) {
      const int m = (l + r) / 2;
      if (nums[m] == target)
        return m;
      if (nums[m] < target)
        l = m + 1;
      else
        r = m;
    }

    return l;
  }
};
class Solution {
 public:
  int lengthOfLastWord(string s) {
    int i = s.length() - 1;

    while (i >= 0 && s[i] == ' ')
      --i;
    const int lastIndex = i;
    while (i >= 0 && s[i] != ' ')
      --i;

    return lastIndex - i;
  }
};

class Solution {
 public:
  vector<int> plusOne(vector<int>& digits) {
    for (int i = digits.size() - 1; i >= 0; --i) {
      if (digits[i] < 9) {
        ++digits[i];
        return digits;
      }
      digits[i] = 0;
    }

    digits.insert(digits.begin(), 1);
    return digits;
  }
};
class Solution {
 public:
  string addBinary(string a, string b) {
    string ans;
    int carry = 0;
    int i = a.length() - 1;
    int j = b.length() - 1;

    while (i >= 0 || j >= 0 || carry) {
      if (i >= 0)
        carry += a[i--] - '0';
      if (j >= 0)
        carry += b[j--] - '0';
      ans += carry % 2 + '0';
      carry /= 2;
    }

    ranges::reverse(ans);
    return ans;
  }
};

class Solution {
 public:
  int mySqrt(int x) {
    unsigned l = 1;
    unsigned r = x + 1u;

    while (l < r) {
      const unsigned m = (l + r) / 2;
      if (m > x / m)
        r = m;
      else
        l = m + 1;
    }

    // l := the minimum number s.t. l * l > x
    return l - 1;
  }
};
class Solution {
 public:
  int climbStairs(int n) {
    // dp[i] := the number of ways to climb to the i-th stair
    vector<int> dp(n + 1);
    dp[0] = 1;
    dp[1] = 1;
    class Solution {
 public:
  void merge(vector<int>& nums1, int m, vector<int>& nums2, int n) {
    int i = m - 1;      // nums1's index (the actual nums)
    int j = n - 1;      // nums2's index
    int k = m + n - 1;  // nums1's index (the next filled position)

    while (j >= 0)
      if (i >= 0 && nums1[i] > nums2[j])
        nums1[k--] = nums1[i--];
      else
        nums1[k--] = nums2[j--];
  }
};
class Solution {
 public:
  void merge(vector<int>& nums1, int m, vector<int>& nums2, int n) {
    int i = m - 1;      // nums1's index (the actual nums)
    int j = n - 1;      // nums2's index
    int k = m + n - 1;  // nums1's index (the next filled position)

    while (j >= 0)
      if (i >= 0 && nums1[i] > nums2[j])
        nums1[k--] = nums1[i--];
      else
        nums1[k--] = nums2[j--];
  }
};
class Solution {
 public:
  bool isSameTree(TreeNode* p, TreeNode* q) {
    if (!p || !q)
      return p == q;
    return p->val == q->val &&              //
           isSameTree(p->left, q->left) &&  //
           isSameTree(p->right, q->right);
  }
};

class Solution {
 public:
  int maxDepth(TreeNode* root) {
    if (root == nullptr)
      return 0;
    return 1 + max(maxDepth(root->left), maxDepth(root->right));
  }
};
lass Solution {
 public:
  TreeNode* sortedArrayToBST(vector<int>& nums) {
    return build(nums, 0, nums.size() - 1);
  }

 private:
  TreeNode* build(const vector<int>& nums, int l, int r) {
    if (l > r)
      return nullptr;
    const int m = (l + r) / 2;
    return new TreeNode(nums[m], build(nums, l, m - 1), build(nums, m + 1, r));
  }
};
lass Solution {
 public:
  vector<vector<int>> generate(int numRows) {
    vector<vector<int>> ans;

    for (int i = 0; i < numRows; ++i)
      ans.push_back(vector<int>(i + 1, 1));

    for (int i = 2; i < numRows; ++i)
      for (int j = 1; j < ans[i].size() - 1; ++j)
        ans[i][j] = ans[i - 1][j - 1] + ans[i - 1][j];

    return ans;
  }
};
class Solution {
 public:
  int maxProfit(vector<int>& prices) {
    int sellOne = 0;
    int holdOne = INT_MIN;

    for (const int price : prices) {
      sellOne = max(sellOne, holdOne + price);
      holdOne = max(holdOne, -price);
    }

    return sellOne;
  }
};

class Solution {
 public:
  bool isPalindrome(string s) {
    int l = 0;
    int r = s.length() - 1;

    while (l < r) {
      while (l < r && !isalnum(s[l]))
        ++l;
      while (l < r && !isalnum(s[r]))
        --r;
      if (tolower(s[l]) != tolower(s[r]))
        return false;
      ++l;
      --r;
    }

    return true;
  }
};
lass Solution {
 public:
  int singleNumber(vector<int>& nums) {
    int ans = 0;

    for (const int num : nums)
      ans ^= num;

    return ans;
  }
};
class Solution {
 public:
  TreeNode* sortedArrayToBST(vector<int>& nums) {
    return build(nums, 0, nums.size() - 1);
  }

 private:
  TreeNode* build(const vector<int>& nums, int l, int r) {
    if (l > r)
      return nullptr;
    const int m = (l + r) / 2;
    return new TreeNode(nums[m], build(nums, l, m - 1), build(nums, m + 1, r));
  }
};
int[] nums = [...]; // Input array
int[] expectedNums = [...]; // The expected answer with correct length

int k = removeDuplicates(nums); // Calls your implementation

assert k == expectedNums.length;
for (int i = 0; i < k; i++) {
    assert nums[i] == expectedNums[i];
}
class Solution {
 public:
  bool hasCycle(ListNode* head) {
    ListNode* slow = head;
    ListNode* fast = head;

    while (fast != nullptr && fast->next != nullptr) {
      slow = slow->next;
      fast = fast->next->next;
      if (slow == fast)
        return true;
    }

    return false;
  }
};
class Solution {
 public:
  vector<int> preorderTraversal(TreeNode* root) {
    vector<int> ans;
    preorder(root, ans);
    return ans;
  }

 private:
  void preorder(TreeNode* root, vector<int>& ans) {
    if (root == nullptr)
      return;

    ans.push_back(root->val);
    preorder(root->left, ans);
    preorder(root->right, ans);
  }
};
class Solution {
 public:
  vector<int> postorderTraversal(TreeNode* root) {
    vector<int> ans;
    postorder(root, ans);
    return ans;
  }

 private:
  void postorder(TreeNode* root, vector<int>& ans) {
    if (root == nullptr)
      return;

    postorder(root->left, ans);
    postorder(root->right, ans);
    ans.push_back(root->val);
  }
};
class Solution {
 public:
  string convertToTitle(int n) {
    return n == 0 ? ""
                  : convertToTitle((n - 1) / 26) + (char)('A' + ((n - 1) % 26));
  }
};
class Solution {
 public:
  ListNode* getIntersectionNode(ListNode* headA, ListNode* headB) {
    ListNode* a = headA;
    ListNode* b = headB;

    while (a != b) {
      a = a == nullptr ? headB : a->next;
      b = b == nullptr ? headA : b->next;
    }

    return a;
  }
};
class Solution {
 public:
  int majorityElement(vector<int>& nums) {
    int ans;
    int count = 0;

    for (const int num : nums) {
      if (count == 0)
        ans = num;
      count += num == ans ? 1 : -1;
    }

    return ans;
  }
};
class Solution {
 public:
  int titleToNumber(string columnTitle) {
    return accumulate(
        columnTitle.begin(), columnTitle.end(), 0,
        [](int subtotal, char c) { return subtotal * 26 + (c - 'A' + 1); });
  }
};
SELECT
  Person.firstName,
  Person.lastName,
  Address.city,
  Address.state
FROM Person
LEFT JOIN Address
  USING (personId);
  SELECT Customers.name AS Customers
FROM Customers
LEFT JOIN Orders
  ON (Customers.id = Orders.customerId)
WHERE Orders.id IS NULL;
