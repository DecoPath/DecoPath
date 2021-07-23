from django.template.loader import render_to_string
from django.utils.encoding import force_bytes
from django.utils.http import urlsafe_base64_encode


from viewer.models import User
from viewer.src.response_handler import *
from viewer.tasks import email
from viewer.tokens import account_activation_token


def verify_email(domain: str, user: User, register_email: str):
    uid = urlsafe_base64_encode(force_bytes(user.email))
    token = account_activation_token.make_token(user)

    act_url = f"http://{domain}/activate/{uid}/{token}"

    html_content = render_to_string(
        'viewer/email_template.html',
        context={"register_email": register_email, "act_url": act_url}
    )
    text_content = render_to_string(
        'viewer/email_template.txt',
        context={"register_email": register_email, "act_url": act_url}
    )

    email.delay(
        subject=EMAIL_SUBJECT,
        text_content=text_content,
        html_content=html_content,
        sender=SENDER_EMAIL,
        recipient_list=[register_email]
    )
    return None
